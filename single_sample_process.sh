#!/bin/bash
# ==============================================================================
# 单样本高并发生信流水线 - v3.1 (全链路降采样极速版)
# 特性：HUMAnN 接入 500 万对降采样流，彻底消灭内存盘 OOM 死锁，极致提速
# ==============================================================================

# ---------- 0. 参数接收与默认配置 ----------
SAMPLE=$1
R1=$2
R2=$3
THREADS=${4:-20}

# 模块开关
RUN_STEP1_CENT=${RUN_STEP1_CENT:-true}
RUN_STEP2_FASTP=${RUN_STEP2_FASTP:-true} 
RUN_STEP3_MPA=${RUN_STEP3_MPA:-true}
RUN_STEP35_STRAIN=${RUN_STEP35_STRAIN:-true}
RUN_STEP4_HUMANN=${RUN_STEP4_HUMANN:-true}
RUN_STEP5_ARGS=${RUN_STEP5_ARGS:-true}
RUN_STEP6_INSTRAIN=${RUN_STEP6_INSTRAIN:-true} 
RUN_STEP7_SGVF=${RUN_STEP7_SGVF:-true}  
RUN_CLEANUP=${RUN_CLEANUP:-true}


# ---------- 1. 路径与环境配置 ----------
source "/data/wangjiaxuan/miniconda3/etc/profile.d/conda.sh"
conda activate meta
export PYTHONWARNINGS="ignore"
set -o pipefail  # 极其重要！

CENT_BIN="/data/wangjiaxuan/miniconda3/envs/meta/bin/centrifuger"
CENT_DB="/data/wangjiaxuan/refer/centrifuger_db/CHM13/chm13_index"
DB_PKL="/data/wangjiaxuan/refer/metaphlan_db/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl"

CLEAN_DIR="1.Clean_Data"
QUANT_DIR="2.Quantity"
ARG_DIR="3.Annotation/ARG"
LOG_DIR="logs"
mkdir -p $CLEAN_DIR $QUANT_DIR $ARG_DIR $LOG_DIR

# (虽然去掉了HUMAnN的内存盘，但为了兼容之前的逻辑，保留清理代码)
SHM_DIR="/dev/shm/hmn_${SAMPLE}"
trap 'rm -rf "$SHM_DIR"' EXIT

LOG_FILE="$LOG_DIR/${SAMPLE}.process.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ---------- 2. 工具函数 ----------
function get_time() {
    date "+%Y-%m-%d %H:%M:%S"
}

function timestamp_start() {
    MODULE_NAME=$1
    CUR_START=$(date +%s)
    echo "[$(get_time)] >>> START: $MODULE_NAME"
}

function timestamp_end() {
    CUR_END=$(date +%s)
    DIFF=$((CUR_END - CUR_START))
    NOW_TIME=$(get_time)
    
    echo "[$NOW_TIME] <<< END: $MODULE_NAME | Time: ${DIFF}s"
    echo "[TIME_AUDIT] [$NOW_TIME] Sample: $SAMPLE | Module: $MODULE_NAME | Time: ${DIFF}s" >> $LOG_DIR/time_summary.audit
}

function check_disk_safety() {
    FREE_MB=$(df -m /data | tail -1 | awk '{print $4}')
    if [[ "$FREE_MB" -lt 153600 ]]; then
        echo "[$(get_time)] Disk space low ($FREE_MB MB), waiting 10m..."
        sleep 600
    fi
    return 0 
}

# ---------- 执行流程 ----------
echo "======= [$(get_time)] Processing: $SAMPLE ======="

if [[ -f "$LOG_DIR/${SAMPLE}.all_done" ]]; then
    echo "[INFO] Sample $SAMPLE is already marked as DONE. Skipping entirely."
    exit 0
fi

# ==============================================================================
# --- Step 1: 粗筛 (fastp 初始质控 + Centrifuger 极速去宿主) ---
# ==============================================================================
QC_IN1="${CLEAN_DIR}/${SAMPLE}.qc_1.fq.gz"
QC_IN2="${CLEAN_DIR}/${SAMPLE}.qc_2.fq.gz"

CENT_UN_PREFIX="${CLEAN_DIR}/${SAMPLE}_cent_unmapped"
CENT_UN1="${CENT_UN_PREFIX}_1.fq.gz"
CENT_UN2="${CENT_UN_PREFIX}_2.fq.gz"
JSON_RAW="${LOG_DIR}/${SAMPLE}.fastp_raw.json"

if [[ "$RUN_STEP1_CENT" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s1.done" ]]; then
    check_disk_safety
    timestamp_start "Step1_Coarse_Filter (fastp + Centrifuger)"

    fastp -i "$R1" -I "$R2" \
          -o "$QC_IN1" -O "$QC_IN2" \
          -q 20 -u 30 -l 50 \
          -w 8 \
          -j "$JSON_RAW" -h /dev/null

    if [[ ! -f "$QC_IN1" ]]; then
        echo "[ERROR] Step 1.1 fastp failed."
        exit 1
    fi

    if $CENT_BIN -x $CENT_DB -1 "$QC_IN1" -2 "$QC_IN2" -t "$THREADS" --min-hitlen 25 \
        --un "$CENT_UN_PREFIX" > /dev/null; then
        
        rm -f "$QC_IN1" "$QC_IN2" 
        touch "$LOG_DIR/${SAMPLE}.s1.done"
        timestamp_end
    else
        echo "[ERROR] Step 1.2 Centrifuger failed."
        exit 1
    fi
fi

# ==============================================================================
# --- Step 2: 精筛 (Bowtie2 深度去宿主 + 二次 fastp + 极速报表) ---
# ==============================================================================
BOWTIE2_DB="/data/wangjiaxuan/refer/kneaddata_db/hg_39"
BT2_UN1="${CLEAN_DIR}/${SAMPLE}.bt2_clean_1.fq.gz"
BT2_UN2="${CLEAN_DIR}/${SAMPLE}.bt2_clean_2.fq.gz"

OUT1="${CLEAN_DIR}/${SAMPLE}.clean_1.fastq"
OUT2="${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"
JSON_CLN="${LOG_DIR}/${SAMPLE}.fastp_clean.json"

if [[ "$RUN_STEP2_FASTP" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s2.done" ]]; then
    timestamp_start "Step2_Fine_Filter (Bowtie2 + fastp_Report)"

    if [[ ! -s "$CENT_UN1" ]]; then
        echo "[WARNING] Centrifuger outputs are completely empty."
        touch "$OUT1" "$OUT2"
        echo '{"summary":{"before_filtering":{"total_reads":0,"q30_rate":0,"gc_content":0},"after_filtering":{"total_reads":0,"q30_rate":0,"gc_content":0}}}' > "$JSON_CLN"
    else
        bowtie2 -p 8 -x "$BOWTIE2_DB" \
                -1 "$CENT_UN1" -2 "$CENT_UN2" \
                --very-sensitive \
                --dovetail \
                --mm \
                --un-conc-gz "${CLEAN_DIR}/${SAMPLE}.bt2_clean_%.fq.gz" \
                > /dev/null 2> "${LOG_DIR}/${SAMPLE}.bowtie2_host.log"
                
        fastp -i "${CLEAN_DIR}/${SAMPLE}.bt2_clean_1.fq.gz" -I "${CLEAN_DIR}/${SAMPLE}.bt2_clean_2.fq.gz" \
              -o "$OUT1" -O "$OUT2" \
              -w 8 \
              -j "$JSON_CLN" -h /dev/null
    fi

    if [[ -f "$JSON_RAW" ]] && [[ -f "$JSON_CLN" ]]; then
        STATS_RES=$(python -c "
import json, sys
try:
    with open('$JSON_RAW') as f1, open('$JSON_CLN') as f2:
        j1 = json.load(f1)
        j2 = json.load(f2)
        r_reads = j1['summary']['before_filtering']['total_reads']
        r_q30 = j1['summary']['before_filtering']['q30_rate'] * 100
        r_gc = j1['summary']['before_filtering']['gc_content'] * 100
        c_reads = j2['summary']['after_filtering']['total_reads']
        c_q30 = j2['summary']['after_filtering']['q30_rate'] * 100
        c_gc = j2['summary']['after_filtering']['gc_content'] * 100
        print(f'{r_reads}\t{r_q30:.2f}%\t{r_gc:.2f}%\t{c_reads}\t{c_q30:.2f}%\t{c_gc:.2f}%')
except Exception as e:
    print('0\t0.00%\t0.00%\t0\t0.00%\t0.00%')
")
        REPORT_FILE="$LOG_DIR/QC_Global_Report.tsv"
        if [[ ! -f "$REPORT_FILE" ]]; then
            echo -e "SampleID\tRaw_Reads\tRaw_Q30\tRaw_GC\tClean_Reads\tClean_Q30\tClean_GC" > "$REPORT_FILE"
        fi
        echo -e "${SAMPLE}\t${STATS_RES}" >> "$REPORT_FILE"
        
        echo "[Cleanup] Removing intermediate Host sequences and JSON files..."
        rm -f "$CENT_UN1" "$CENT_UN2" "${CLEAN_DIR}/${SAMPLE}.bt2_clean_1.fq.gz" "${CLEAN_DIR}/${SAMPLE}.bt2_clean_2.fq.gz" "$JSON_RAW" "$JSON_CLN"
        
        echo "[INFO] QC & Host Depletion finished. Preparing merged Fastq for downstream..."
        CLEAN_FQ_COMBINED="${CLEAN_DIR}/${SAMPLE}.clean.fastq"
        cat "$OUT1" "$OUT2" > "$CLEAN_FQ_COMBINED"
        
        touch "$LOG_DIR/${SAMPLE}.s2.done"
        timestamp_end
    else
        echo "[ERROR] Missing JSON files for QC parsing."
        exit 1
    fi
fi

# ==============================================================================
# --- 核心优化：数据状态自检与“按需”恢复 ---
# ==============================================================================
if [[ -f "$LOG_DIR/${SAMPLE}.s2.done" ]]; then
    NEED_PAIRED="false" 
    NEED_MERGED="false" 

    [[ "$RUN_STEP3_MPA" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s3.done" ]] && NEED_MERGED="true"
    [[ "$RUN_STEP4_HUMANN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]] && NEED_PAIRED="true" # 现在HUMAnN需要双端去合并
    [[ "$RUN_STEP5_ARGS" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]] && NEED_PAIRED="true"
    [[ "$RUN_STEP6_INSTRAIN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s6.done" ]] && NEED_PAIRED="true"
    [[ "$RUN_STEP7_SGVF" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s7.done" ]] && NEED_PAIRED="true" 
    
    [[ "$NEED_MERGED" == "true" ]] && NEED_PAIRED="true"

    if [[ "$NEED_PAIRED" == "true" ]]; then
        for suffix in "1" "2"; do
            FQ_GZ="${CLEAN_DIR}/${SAMPLE}.clean_${suffix}.fastq.gz"
            FQ_UNCOMP="${CLEAN_DIR}/${SAMPLE}.clean_${suffix}.fastq"
            if [[ -f "$FQ_GZ" && ! -f "$FQ_UNCOMP" ]]; then
                echo "[Recovery] Downstream tasks require $FQ_UNCOMP. Decompressing..."
                pigz -d -p "$THREADS" "$FQ_GZ"
            fi
        done
    fi

    CLEAN_FQ_COMBINED="${CLEAN_DIR}/${SAMPLE}.clean.fastq"
    if [[ "$NEED_MERGED" == "true" && ! -f "$CLEAN_FQ_COMBINED" && -f "${CLEAN_DIR}/${SAMPLE}.clean_1.fastq" ]]; then
        echo "[Recovery] Downstream tasks require merged fastq. Re-merging..."
        cat "${CLEAN_DIR}/${SAMPLE}.clean_1.fastq" "${CLEAN_DIR}/${SAMPLE}.clean_2.fastq" > "$CLEAN_FQ_COMBINED"
    fi
fi

CLEAN_FQ="${CLEAN_DIR}/${SAMPLE}.clean.fastq"

# --- Step 3: MetaPhlAn 4 (带失败自动清理) ---
# 注意：MetaPhlAn 依然使用抽样前的全量数据 $CLEAN_FQ
if [[ "$RUN_STEP3_MPA" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s3.done" ]]; then
    timestamp_start "MetaPhlAn"
    metaphlan "$CLEAN_FQ" -o "${QUANT_DIR}/${SAMPLE}_profile.txt" \
              --input_type fastq --nproc "$THREADS" --offline --sample_id "${SAMPLE}" \
              -x mpa_vOct22_CHOCOPhlAnSGB_202403 -s "${QUANT_DIR}/${SAMPLE}.sam.bz2" \
              --mapout "${QUANT_DIR}/${SAMPLE}.bt2out" \
              -t rel_ab_w_read_stats \
              --profile_vsc \
              --vsc_out "${QUANT_DIR}/${SAMPLE}_vsc_profile.txt"

    if [[ -s "${QUANT_DIR}/${SAMPLE}_profile.txt" ]]; then
        rm -f "${QUANT_DIR}/${SAMPLE}.bt2out" 
        touch "$LOG_DIR/${SAMPLE}.s3.done" ; timestamp_end
    else
        echo "[ERROR] MetaPhlAn failed. Cleaning up residuals..."
        rm -f "${QUANT_DIR}/${SAMPLE}_profile.txt" "${QUANT_DIR}/${SAMPLE}.bt2out" "${QUANT_DIR}/${SAMPLE}.sam.bz2"
        exit 1 
    fi
fi

# --- Step 3.5: StrainPhlAn (Markers) ---
if [[ "$RUN_STEP35_STRAIN" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s3_5.done" ]]; then
    timestamp_start "Strain_Markers"
    mkdir -p "${QUANT_DIR}/db_markers"
    
    if [[ -f "${QUANT_DIR}/db_markers/${SAMPLE}.json.bz2" ]]; then
        echo "[Strain_Markers] Output .json.bz2 already exists. Recovering done mark..."
        [[ -f "${QUANT_DIR}/${SAMPLE}.sam.bz2" ]] && rm -f "${QUANT_DIR}/${SAMPLE}.sam.bz2"
        touch "$LOG_DIR/${SAMPLE}.s3_5.done"
        timestamp_end
    elif [[ -f "${QUANT_DIR}/${SAMPLE}.sam.bz2" ]]; then
        sample2markers.py -i "${QUANT_DIR}/${SAMPLE}.sam.bz2" -o "${QUANT_DIR}/db_markers/" \
                          -n "$THREADS" -f bz2 -d "$DB_PKL" -b 80 \
                          --min_reads_aligning 8
        
        if [[ -f "${QUANT_DIR}/db_markers/${SAMPLE}.json.bz2" ]]; then
            rm -f "${QUANT_DIR}/${SAMPLE}.sam.bz2"
            touch "$LOG_DIR/${SAMPLE}.s3_5.done"
            timestamp_end
        else
            echo "[ERROR] StrainPhlAn marker extraction failed or output format changed."
            exit 1
        fi
    else
        echo "[ERROR] Missing both ${SAMPLE}.sam.bz2 and ${SAMPLE}.json.bz2. Cannot run Step 3.5!"
        exit 1
    fi
fi


# ==============================================================================
# --- 【核心上移】Step 3.8: 全局统一降采样 (智能按需+零延时优化版) ---
# 保护对象：HUMAnN, ARGs-OAP, inStrain, SGVFinder2
# ==============================================================================
MAP_IN1="${CLEAN_DIR}/${SAMPLE}.clean_1.fastq"
MAP_IN2="${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"
SUB_R1="${CLEAN_DIR}/${SAMPLE}.sub_1.fastq"
SUB_R2="${CLEAN_DIR}/${SAMPLE}.sub_2.fastq"

# 设定安全红线：600 万对。抽样目标：500 万对
MAX_PAIRS=6000000
TARGET_PAIRS=5000000
SEQTK_BIN="/data/wangjiaxuan/miniconda3/bin/seqtk"
JSON_CLN="${LOG_DIR}/${SAMPLE}.fastp_clean.json"

# 【神级修复】：智能预判！不仅要看开关是否打开，还要看它是不是还没跑完（缺少 .done）！
NEED_SUBSAMPLE="false"
[[ "$RUN_STEP4_HUMANN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]] && NEED_SUBSAMPLE="true"
[[ "$RUN_STEP5_ARGS" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]] && NEED_SUBSAMPLE="true"
[[ "$RUN_STEP6_INSTRAIN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s6.done" ]] && NEED_SUBSAMPLE="true"
[[ "$RUN_STEP7_SGVF" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s7.done" ]] && NEED_SUBSAMPLE="true" 

# 只有当下游真真切切需要跑的时候，才进入抽样逻辑
if [[ "$NEED_SUBSAMPLE" == "true" ]]; then
    
    # 优化 1：如果曾经抽样过且文件还在，直接复用，跳过一切计算！
    if [[ -s "$SUB_R1" && -s "$SUB_R2" ]]; then
        echo "[Pre-Processing] Subsampled files already exist. Reusing them instantly..."
        MAP_IN1="$SUB_R1"
        MAP_IN2="$SUB_R2"
    else
        # 优化 2：从 fastp 的 JSON 里秒读 reads 数，砍掉 wc -l 的 5 分钟卡顿！
        if [[ -f "$JSON_CLN" ]]; then
            CUR_PAIRS=$(grep -o '"total_reads":[^,]*' "$JSON_CLN" | head -1 | awk -F':' '{print int($2/2)}')
        else
            CUR_PAIRS=$(($(wc -l < "$MAP_IN1") / 4))
        fi
        
        if [ "$CUR_PAIRS" -gt "$MAX_PAIRS" ]; then
            echo "[Pre-Processing] Sample depth is high (${CUR_PAIRS} pairs). Subsampling to ${TARGET_PAIRS}..."
            
            if ! $SEQTK_BIN sample -s 100 "$MAP_IN1" "$TARGET_PAIRS" > "${SUB_R1}.tmp"; then
                echo "[ERROR] seqtk failed on R1!" && exit 1
            fi
            if ! $SEQTK_BIN sample -s 100 "$MAP_IN2" "$TARGET_PAIRS" > "${SUB_R2}.tmp"; then
                echo "[ERROR] seqtk failed on R2!" && exit 1
            fi
            mv "${SUB_R1}.tmp" "$SUB_R1"
            mv "${SUB_R2}.tmp" "$SUB_R2"
            
            # 更新指针
            MAP_IN1="$SUB_R1"
            MAP_IN2="$SUB_R2"
        else
            echo "[Pre-Processing] Sample depth (${CUR_PAIRS} pairs) is normal. No subsampling needed."
        fi
    fi
else
    # 【新增逻辑反馈】：如果下游都跑完了或者没开，清楚地告诉用户跳过抽样
    echo "[Pre-Processing] All downstream modules are completed or disabled. Skipping subsampling."
fi

# ==============================================================================
# --- Step 4: HUMAnN 3 (终极形态：降采样输入 + 内存盘极速 + 动态排队红绿灯) ---
# ==============================================================================
if [[ "$RUN_STEP4_HUMANN" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]]; then
    timestamp_start "HUMAnN"
    
    # 1. 专门为 HUMAnN 合并抽样后的 Clean Reads
    SUB_MERGED="${CLEAN_DIR}/${SAMPLE}.sub_merged.fastq"
    if [ ! -s "$SUB_MERGED" ]; then
        echo "[HUMAnN] Merging subsampled reads for HUMAnN input..."
        cat "$MAP_IN1" "$MAP_IN2" > "$SUB_MERGED"
    fi

    # ==============================================================================
    # 2. 内存盘动态“红绿灯”排队机制
    # ==============================================================================
    # 因为输入固定为 500 万对，Diamond 极其轻量，这里门槛可以放心降到 20GB！
    REQUIRED_GB=20  
    
    # 随机错峰休眠 (1~60秒)，防止多个样本同秒惊群抢占
    sleep $((RANDOM % 60)) 
    
    while true; do
        FREE_SHM=$(df -BG /dev/shm | tail -1 | awk '{print $4}' | sed 's/G//')
        
        if [[ "$FREE_SHM" -ge "$REQUIRED_GB" ]]; then
            echo "[$(date)] [HUMAnN] /dev/shm has ${FREE_SHM}GB free. Space acquired for $SAMPLE."
            break # 空间充足，跳出排队，开始干活！
        else
            echo "[$(date)] [HUMAnN] /dev/shm full (Free: ${FREE_SHM}GB < ${REQUIRED_GB}GB). $SAMPLE waiting 5 minutes..."
            sleep 300 # 空间不足，睡 5 分钟后再检查
        fi
    done
    # ==============================================================================

    # 3. 杀入内存盘极速运行 (自带 timeout 6h 终极防死锁保险)
    timeout 6h humann --threads "$THREADS" \
           --input "$SUB_MERGED" \
           --output "$SHM_DIR" \
           --output-basename "$SAMPLE" \
           --taxonomic-profile "${QUANT_DIR}/${SAMPLE}_profile.txt" \
           --remove-temp-output
           
    EXIT_CODE=$?
    
    # 4. 把跑完的丰度表从内存盘拷回物理硬盘
    cp ${SHM_DIR}/${SAMPLE}_*.[tl]* "$QUANT_DIR/" 2>/dev/null

    # 5. 验证是否成功与极限清理
    if [[ $EXIT_CODE -eq 124 ]]; then
        echo "[ERROR] HUMAnN timeout (>6h), possible lock. Killed by system."
        rm -rf "${SHM_DIR:?}"/*
        rm -f "$SUB_MERGED"
        exit 1
    elif [[ -f "${QUANT_DIR}/${SAMPLE}_genefamilies.tsv" ]] || [[ -f "${QUANT_DIR}/${SAMPLE}_2_genefamilies.tsv" ]]; then
        touch "$LOG_DIR/${SAMPLE}.s4.done"
        
        # 【关键清理】：立刻清空当前样本的内存盘，马上把内存还给服务器！合并的中间文件也删掉！
        rm -rf "${SHM_DIR:?}"/*
        rm -f "$SUB_MERGED" 
        
        timestamp_end
    else
        echo "[ERROR] HUMAnN failed with exit code $EXIT_CODE. Check logs."
        rm -rf "${SHM_DIR:?}"/*
        exit 1
    fi
fi

# --- Step 5: ARGs-OAP (软链接共享抽样版) ---
if [[ "$RUN_STEP5_ARGS" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]]; then
    timestamp_start "ARGs_OAP"
    W_DIR="${ARG_DIR}/${SAMPLE}" ; mkdir -p "${W_DIR}/work_space"
    
    echo "[ARGs-OAP] Linking global fastq files to work_space (Zero I/O overhead)..."
    ln -sf "$(realpath "$MAP_IN1")" "${W_DIR}/work_space/sample_1.fastq"
    ln -sf "$(realpath "$MAP_IN2")" "${W_DIR}/work_space/sample_2.fastq"

    PUSH_DIR=$(pwd)
    cd "${W_DIR}/work_space" || exit

    args_oap stage_one -i ./ -o ./result -f fastq -t "$THREADS"
    args_oap stage_two -i ./result -t "$THREADS" -o ./result
    
    if [[ -d "./result" ]]; then
        mv ./result/normalized_cell.*.txt ../
        cd "$PUSH_DIR" || exit
        rm -rf "${W_DIR}/work_space" 
        [[ -f "${PUSH_DIR}/fifo_map.mapout.txt" ]] && rm -f "${PUSH_DIR}/fifo_map.mapout.txt"
        touch "$LOG_DIR/${SAMPLE}.s5.done" ; timestamp_end
    else
        cd "$PUSH_DIR" || exit
        echo "[ERROR] ARGs-OAP failed. Wiping work_space..."
        rm -rf "${W_DIR}/work_space" 
        [[ -f "${PUSH_DIR}/fifo_map.mapout.txt" ]] && rm -f "${PUSH_DIR}/fifo_map.mapout.txt"
        exit 1
    fi
fi

# --- Step 6: inStrain (微多样性与群体遗传学全景分析) ---
IS_DB_DIR="/data/wangjiaxuan/refer/preterm_refseq_400/preterm_refseq_400_v1"
IS_FASTA="${IS_DB_DIR}/preterm_refseq_400.fasta"
IS_BT2="${IS_DB_DIR}/preterm_refseq_400.fasta.bt2"
IS_STB="${IS_DB_DIR}/preterm_refseq_400.stb"
IS_GENES="${IS_DB_DIR}/preterm_refseq_400.genes.fna"

if [[ "$RUN_STEP6_INSTRAIN" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s6.done" ]]; then
    timestamp_start "inStrain"
    
    BAM_OUT="${QUANT_DIR}/${SAMPLE}_preterm400.sorted.bam"
    IS_OUT_DIR="${QUANT_DIR}/${SAMPLE}_inStrain.IS"

    if [[ ! -f "$BAM_OUT" ]]; then
        echo "[inStrain] Running Bowtie2 mapping to preterm_refseq_400 database..."
        bowtie2 -p "$THREADS" -x "$IS_BT2" -1 "$MAP_IN1" -2 "$MAP_IN2" --no-unal | \
        samtools sort -@ "$THREADS" -m 2G -o "$BAM_OUT"
        
        echo "[inStrain] Indexing BAM file..."
        samtools index -@ "$THREADS" "$BAM_OUT"
    fi

    if [[ -f "$BAM_OUT" ]]; then
        echo "[inStrain] Running inStrain profile with full parameters..."
        /data/wangjiaxuan/miniconda3/envs/meta1/bin/inStrain  \
        profile "$BAM_OUT" "$IS_FASTA" \
            -o "$IS_OUT_DIR" \
            -p "$THREADS" \
            -g "$IS_GENES" \
            -s "$IS_STB" \
            --database_mode
        
        if [[ -d "${IS_OUT_DIR}/output" ]]; then
            touch "$LOG_DIR/${SAMPLE}.s6.done"
            timestamp_end
            
            if [[ -d "${QUANT_DIR}/${SAMPLE}_inStrain.IS" ]]; then
                echo "[inStrain] Cleaning up bloated logs and figures..."
                rm -rf "${QUANT_DIR}/${SAMPLE}_inStrain.IS/figures"
                rm -rf "${QUANT_DIR}/${SAMPLE}_inStrain.IS/log"
                rm -f "$BAM_OUT" "${BAM_OUT}.bai"
            fi
        else
            echo "[ERROR] inStrain profile failed to generate output directory. Check logs."
            exit 1
        fi
    else
        echo "[ERROR] BAM file was not generated successfully. Aborting inStrain."
        exit 1
    fi
fi

# --- Step 7: SGVFinder2 (大片段结构变异与基因缺失分析) ---
if [[ "$RUN_STEP7_SGVF" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s7.done" ]]; then
    timestamp_start "SGVFinder2"
    
    ICRA_BIN="/data/wangjiaxuan/miniconda3/envs/sgvfinder/bin/icra"
    SVF_BIN="/data/wangjiaxuan/miniconda3/envs/sgvfinder/bin/svfinder"
    SGVF_DB_PREFIX="/data/wangjiaxuan/refer/preterm_refseq_400/sgvfinder/preterm_refseq_400_sgvfinder"
    SGVF_OUT_DIR="${QUANT_DIR}/${SAMPLE}_SGVF"
    
    mkdir -p "$SGVF_OUT_DIR"
    
    JSDEL_FILE=$(find "$SGVF_OUT_DIR" -maxdepth 1 -name "*.jsdel" | head -n 1)

    if [[ -z "$JSDEL_FILE" ]]; then
        echo "[SGVFinder2] Running ICRA (Iterative Coverage-based Read Assignment)..."
        rm -f "$SGVF_OUT_DIR"/*.bam
        
        $ICRA_BIN \
            --fq1 "$MAP_IN1" \
            --fq2 "$MAP_IN2" \
            --db "$SGVF_DB_PREFIX" \
            --outfol "$SGVF_OUT_DIR" \
            --threads "$THREADS"
            
        JSDEL_FILE=$(find "$SGVF_OUT_DIR" -maxdepth 1 -name "*.jsdel" | head -n 1)
    fi

    if [[ -n "$JSDEL_FILE" ]]; then
        echo "[SGVFinder2] Generating sample coverage map from: $(basename "$JSDEL_FILE")"
        
        $SVF_BIN get_sample_map --delta_file "$JSDEL_FILE" --db_path "$SGVF_DB_PREFIX"
        
        if [[ $? -eq 0 ]]; then
            touch "$LOG_DIR/${SAMPLE}.s7.done"
            rm -f "$SGVF_OUT_DIR"/*.bam "$SGVF_OUT_DIR"/*.pmp "$SGVF_OUT_DIR"/*.jspi "$SGVF_OUT_DIR"/*.jsdel
            timestamp_end
        else
            echo "[ERROR] SGVFinder2 get_sample_map failed. Check logs."
            exit 1
        fi
    else
        echo "[ERROR] SGVFinder2 (ICRA) failed to generate .jsdel file. Check logs."
        exit 1
    fi
fi

# --- Final Cleanup (动态感知全通过版) ---
if [[ "$RUN_CLEANUP" == "true" ]]; then
    CAN_CLEAN="true"
    [[ "$RUN_STEP3_MPA" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s3.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP35_STRAIN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s3_5.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP4_HUMANN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP5_ARGS" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP6_INSTRAIN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s6.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP7_SGVF" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s7.done" ]] && CAN_CLEAN="false" 

    if [[ "$CAN_CLEAN" == "true" ]]; then
        echo "[$(get_time)] >>> ALL TARGET STEPS SUCCESSFUL: $SAMPLE. Starting cleanup."
        [[ -f "${QUANT_DIR}/${SAMPLE}.sam.bz2" ]] && rm -f "${QUANT_DIR}/${SAMPLE}.sam.bz2"

        if [[ -f "${CLEAN_DIR}/${SAMPLE}.clean.fastq" ]]; then
            rm -f "${CLEAN_DIR}/${SAMPLE}.clean.fastq"
        fi

        # 彻底销毁统一降采样产生的临时文件
        SUB_R1="${CLEAN_DIR}/${SAMPLE}.sub_1.fastq"
        SUB_R2="${CLEAN_DIR}/${SAMPLE}.sub_2.fastq"
        SUB_MERGED="${CLEAN_DIR}/${SAMPLE}.sub_merged.fastq"
        
        echo "[Cleanup] Removing subsampled temporary fastq files..."
        rm -f "$SUB_R1" "$SUB_R2" "$SUB_MERGED"

        # 压缩双端 clean data 留档
        for fq in "${CLEAN_DIR}/${SAMPLE}.clean_1.fastq" "${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"; do
            if [[ -f "$fq" ]]; then
                echo "Compressing $fq with pigz..."
                pigz -f -p 8 "$fq"
            fi
        done
        
        touch "$LOG_DIR/${SAMPLE}.all_done"
    else
        echo "[$(get_time)] >>> SKIP CLEANUP: $SAMPLE has incomplete steps (or failed). Check logs."
    fi
fi

echo "======= [$(get_time)] FINISHED: $SAMPLE ======="