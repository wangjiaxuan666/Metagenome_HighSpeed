#!/bin/bash
# ==============================================================================
# 单样本高并发生信流水线 - v3.0 (fastp 加速版)
# 特性：用 fastp 替代 Kneaddata，解决 Java 内存溢出问题，大幅提升速度
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
RUN_STEP6_INSTRAIN=${RUN_STEP6_INSTRAIN:-true} # <-- 移到这里统一管理
RUN_STEP7_SGVF=${RUN_STEP7_SGVF:-true}  # <-- 新增这一行
RUN_CLEANUP=${RUN_CLEANUP:-true}


# ---------- 1. 路径与环境配置 ----------
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate meta
export PYTHONWARNINGS="ignore"

CENT_BIN="/data/wangjiaxuan/miniconda3/envs/meta/bin/centrifuger"
CENT_DB="/data/wangjiaxuan/refer/centrifuger_db/CHM13/chm13_index"
DB_PKL="/data/wangjiaxuan/refer/metaphlan_db/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl"

CLEAN_DIR="1.Clean_Data"
QUANT_DIR="2.Quantity"
ARG_DIR="3.Annotation/ARG"
LOG_DIR="logs"
mkdir -p $CLEAN_DIR $QUANT_DIR $ARG_DIR $LOG_DIR

# 定义内存盘路径（每个样本独立文件夹）
SHM_DIR="/dev/shm/hmn_${SAMPLE}"
mkdir -p "$SHM_DIR"
# 确保脚本退出（无论成功失败）都会清理内存盘
trap 'rm -rf "$SHM_DIR"' EXIT


LOG_FILE="$LOG_DIR/${SAMPLE}.process.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ---------- 2. 工具函数 ----------

# 定义统一的时间格式，保持和 ARGs-OAP 等 Python 软件的日志格式完全一致
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
    
    # 打印到当前样本的专属日志
    echo "[$NOW_TIME] <<< END: $MODULE_NAME | Time: ${DIFF}s"
    
    # 【修复重点】：把具体的结束时间 ($NOW_TIME) 写入审计文件
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

# 执行流程开头的提示
echo "======= [$(get_time)] Processing: $SAMPLE ======="

# --- Step 1: Centrifuger (去宿主) ---
# 注意：Centrifuger 出来的文件是 .fq.gz 格式
if [[ "$RUN_STEP1_CENT" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s1.done" ]]; then
    check_disk_safety
    timestamp_start "Centrifuger"
    if $CENT_BIN -x $CENT_DB -1 "$R1" -2 "$R2" -t "$THREADS" --min-hitlen 25 --un "${CLEAN_DIR}/${SAMPLE}_nonhuman" > /dev/null; then
        touch "$LOG_DIR/${SAMPLE}.s1.done"
        timestamp_end
    fi
fi

# --- Step 2: fastp (双端输出 + 兼容合并) ---
if [[ "$RUN_STEP2_FASTP" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s2.done" ]]; then
    timestamp_start "fastp"
    IN1="${CLEAN_DIR}/${SAMPLE}_nonhuman_1.fq.gz"
    IN2="${CLEAN_DIR}/${SAMPLE}_nonhuman_2.fq.gz"
    
    # 定义双端输出，供 inStrain 使用
    OUT1="${CLEAN_DIR}/${SAMPLE}.clean_1.fastq"
    OUT2="${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"
    
    # 取消 --merge，保留标准双端输出
    # 一条 read 修剪后长度不足 50bp，或者它包含过多（>30%）错误率大于 1% 的低质量碱基，就直接把它扔掉。
    fastp -i "$IN1" -I "$IN2" \
          -o "$OUT1" -O "$OUT2" \
          -q 20 -u 30 -l 50 \
          -w "$THREADS" \
          -j "${LOG_DIR}/${SAMPLE}.fastp.json" \
          -h /dev/null

    if [[ -f "$OUT1" ]] && [[ -f "$OUT2" ]]; then
        rm -f "$IN1" "$IN2" # 清理去宿主后的临时文件
        
        # 【关键兼容】：将双端合并为单文件，供 MetaPhlAn / HUMAnN / ARGs-OAP 使用
        CLEAN_FQ_COMBINED="${CLEAN_DIR}/${SAMPLE}.clean.fastq"
        cat "$OUT1" "$OUT2" > "$CLEAN_FQ_COMBINED"
        
        # 将合并后的文件放入内存盘提速
        # cp "$CLEAN_FQ_COMBINED" "$SHM_DIR/${SAMPLE}.clean.fastq"
        # 太占用内存，去掉！

        touch "$LOG_DIR/${SAMPLE}.s2.done"
        timestamp_end
    fi
fi

# ==============================================================================
# --- 核心优化：数据状态自检与“按需”恢复 ---
# ==============================================================================
if [[ -f "$LOG_DIR/${SAMPLE}.s2.done" ]]; then
    # 1. 智能预判：下游到底有没有模块需要跑？需要什么格式的文件？
    NEED_PAIRED="false"  # inStrain 需要双端文件
    NEED_MERGED="false"  # MPA / HUMAnN / ARGs 需要合并版单文件

    [[ "$RUN_STEP3_MPA" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s3.done" ]] && NEED_MERGED="true"
    [[ "$RUN_STEP4_HUMANN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]] && NEED_MERGED="true"
    [[ "$RUN_STEP5_ARGS" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]] && NEED_MERGED="true"
    [[ "$RUN_STEP6_INSTRAIN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s6.done" ]] && NEED_PAIRED="true"
    [[ "$RUN_STEP7_SGVF" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s7.done" ]] && NEED_PAIRED="true" # <-- 新增这一行
    # 逻辑联动：如果需要合并文件，那必然需要先解压双端文件
    [[ "$NEED_MERGED" == "true" ]] && NEED_PAIRED="true"

    # 2. 只有在确定需要双端文件时，才去检查并解压
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

    # 3. 只有在确定需要合并文件时，才去检查并合并
    CLEAN_FQ_COMBINED="${CLEAN_DIR}/${SAMPLE}.clean.fastq"
    if [[ "$NEED_MERGED" == "true" && ! -f "$CLEAN_FQ_COMBINED" && -f "${CLEAN_DIR}/${SAMPLE}.clean_1.fastq" ]]; then
        echo "[Recovery] Downstream tasks require merged fastq. Re-merging..."
        cat "${CLEAN_DIR}/${SAMPLE}.clean_1.fastq" "${CLEAN_DIR}/${SAMPLE}.clean_2.fastq" > "$CLEAN_FQ_COMBINED"
    fi

    # # 4. 只有在确定需要合并文件时，才拷入内存盘
    # if [[ "$NEED_MERGED" == "true" && -f "$CLEAN_FQ_COMBINED" && ! -f "$SHM_DIR/${SAMPLE}.clean.fastq" ]]; then
    #     echo "[Recovery] Copying merged fastq to RAM disk..."
    #     cp "$CLEAN_FQ_COMBINED" "$SHM_DIR/${SAMPLE}.clean.fastq"
    # fi
    # 太占用内存去掉
fi
# ==============================================================================

# 统一规矩：所有输入文件全都在物理硬盘上读，绝不占用宝贵的内存盘！
CLEAN_FQ="${CLEAN_DIR}/${SAMPLE}.clean.fastq"

# --- Step 3: MetaPhlAn 4 (带失败自动清理) ---
if [[ "$RUN_STEP3_MPA" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s3.done" ]]; then
    timestamp_start "MetaPhlAn"
    # 强制落地 bt2out 解决管道 Bug
    metaphlan "$CLEAN_FQ" -o "${QUANT_DIR}/${SAMPLE}_profile.txt" \
              --input_type fastq --nproc "$THREADS" --offline --sample_id "${SAMPLE}" \
              -x mpa_vOct22_CHOCOPhlAnSGB_202403 -s "${QUANT_DIR}/${SAMPLE}.sam.bz2" \
              --mapout "${QUANT_DIR}/${SAMPLE}.bt2out" \
              -t rel_ab_w_read_stats \
              --profile_vsc \
              --vsc_out "${QUANT_DIR}/${SAMPLE}_vsc_profile.txt"

    if [[ -s "${QUANT_DIR}/${SAMPLE}_profile.txt" ]]; then
        rm -f "${QUANT_DIR}/${SAMPLE}.bt2out" # 成功则清理中间文件
        touch "$LOG_DIR/${SAMPLE}.s3.done" ; timestamp_end
    else
        echo "[ERROR] MetaPhlAn failed. Cleaning up residuals..."
        rm -f "${QUANT_DIR}/${SAMPLE}_profile.txt" "${QUANT_DIR}/${SAMPLE}.bt2out" "${QUANT_DIR}/${SAMPLE}.sam.bz2"
        exit 1 # 立即退出，防止后续模块在错误数据上空转
    fi
fi

# --- Step 3.5: StrainPhlAn (Markers) ---
if [[ "$RUN_STEP35_STRAIN" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s3_5.done" ]]; then
    timestamp_start "Strain_Markers"
    mkdir -p "${QUANT_DIR}/db_markers"
    
    # 优化 1：先检查目标结果文件是否已经存在（解决断点死锁问题）
    if [[ -f "${QUANT_DIR}/db_markers/${SAMPLE}.json.bz2" ]]; then
        echo "[Strain_Markers] Output .json.bz2 already exists. Recovering done mark..."
        # 顺手确保多余的 sam 文件被清理
        [[ -f "${QUANT_DIR}/${SAMPLE}.sam.bz2" ]] && rm -f "${QUANT_DIR}/${SAMPLE}.sam.bz2"
        touch "$LOG_DIR/${SAMPLE}.s3_5.done"
        timestamp_end
        
    # 优化 2：如果结果不在，且原料 (sam) 还在，则开始正经干活
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
        
    # 优化 3：既没有结果，也没有原料，直接报错拦截，防止静默失败
    else
        echo "[ERROR] Missing both ${SAMPLE}.sam.bz2 and ${SAMPLE}.json.bz2. Cannot run Step 3.5!"
        exit 1
    fi
fi

# --- Step 4: HUMAnN (带动态排队与内存盘限流保护) ---
if [[ "$RUN_STEP4_HUMANN" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]]; then
    timestamp_start "HUMAnN"
    
    # ==============================================================================
    # 核心优化：内存盘动态“红绿灯”排队机制
    # ==============================================================================
    REQUIRED_GB=50  # 假设每个 HUMAnN 进程至少需要 50GB 空闲内存盘 (可根据你的实际情况调整)
    
    # 1. 随机错峰休眠 (1~60秒)，防止所有刚跑完 MetaPhlAn 的样本在同一秒醒来抢占资源 (Thundering herd problem)
    sleep $((RANDOM % 60)) 
    
    while true; do
        # 2. 实时获取 /dev/shm 的可用空间 (以 GB 为单位)
        FREE_SHM=$(df -BG /dev/shm | tail -1 | awk '{print $4}' | sed 's/G//')
        
        if [[ "$FREE_SHM" -ge "$REQUIRED_GB" ]]; then
            echo "[$(date)] [HUMAnN] /dev/shm has ${FREE_SHM}GB free. Space acquired for $SAMPLE."
            break # 空间充足，跳出循环，开始干活！
        else
            echo "[$(date)] [HUMAnN] /dev/shm full (Free: ${FREE_SHM}GB < ${REQUIRED_GB}GB). $SAMPLE waiting 5 minutes..."
            sleep 300 # 空间不足，睡 5 分钟后再来检查
        fi
    done
    # ==============================================================================

    # 开始运行 HUMAnN
    humann --threads "$THREADS" \
           --input "$CLEAN_FQ" \
           --output "$SHM_DIR" \
           --output-basename "$SAMPLE" \
           --taxonomic-profile "${QUANT_DIR}/${SAMPLE}_profile.txt" \
           --remove-temp-output
           
    # 把结果拷回物理硬盘（完美兼容带有 _2_ _3_ _4_ 前缀的 tsv 和 log）
    cp ${SHM_DIR}/${SAMPLE}_*.[tl]* "$QUANT_DIR/" 2>/dev/null

    # 验证是否成功
    if [[ -f "${QUANT_DIR}/${SAMPLE}_genefamilies.tsv" ]] || [[ -f "${QUANT_DIR}/${SAMPLE}_2_genefamilies.tsv" ]]; then
        touch "$LOG_DIR/${SAMPLE}.s4.done"
        timestamp_end
        
        # 【关键清理】：拷出结果后，立刻清空当前样本的内存盘残留，马上把内存还给服务器！
        rm -rf "${SHM_DIR:?}"/*
    else
        echo "[ERROR] HUMAnN failed. Check logs."
        exit 1
    fi
fi

# --- Step 5: ARGs-OAP (带强制清理逻辑) ---
# 在脚本开头添加抽样参数（比如这里设为单端 5M，双端加起来就是 10M reads）

# --- Step 5: ARGs-OAP (双端同步抽样 + 目录直接读取版) ---
if [[ "$RUN_STEP5_ARGS" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]]; then
    timestamp_start "ARGs_OAP"
    W_DIR="${ARG_DIR}/${SAMPLE}" ; mkdir -p "${W_DIR}/work_space"
    
    # 动态获取 R1 和 R2 的路径 (依赖前面解压/生成的双端文件)
    FQ1="${CLEAN_DIR}/${SAMPLE}.clean_1.fastq"
    FQ2="${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"
    
    SEQTK_BIN="/data/wangjiaxuan/miniconda3/bin/seqtk" 

    echo "[ARGs-OAP] Subsampling to max 5000000 reads per end using seqtk..."
    
    # 直接使用固定值 5000000 抽样
    $SEQTK_BIN sample -s 42 "$FQ1" 5000000 > "${W_DIR}/work_space/sample_1.fastq"
    $SEQTK_BIN sample -s 42 "$FQ2" 5000000 > "${W_DIR}/work_space/sample_2.fastq"

    PUSH_DIR=$(pwd)
    cd "${W_DIR}/work_space" || exit

    # ARGs-OAP 会自动读取当前目录 (./) 下的 sample_1.fastq 和 sample_2.fastq
    args_oap stage_one -i ./ -o ./result -f fastq -t "$THREADS"
    args_oap stage_two -i ./result -t "$THREADS" -o ./result
    
    if [[ -d "./result" ]]; then
        mv ./result/normalized_cell.*.txt ../
        cd "$PUSH_DIR" || exit
        rm -rf "${W_DIR}/work_space" # 成功后整个目录连带抽样的 fastq 一起删掉
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
# 前提：Step 2 fastp 已经修改为输出标准双端 reads ($OUT1 和 $OUT2)

# 1. inStrain 数据库绝对路径定义 (修改为新的 preterm_refseq_400 数据库)
IS_DB_DIR="/data/wangjiaxuan/refer/preterm_refseq_400/preterm_refseq_400_v1"
IS_FASTA="${IS_DB_DIR}/preterm_refseq_400.fasta"
IS_BT2="${IS_DB_DIR}/preterm_refseq_400.fasta.bt2" # Bowtie2 index 前缀 (完美兼容 .bt2l 大索引)
IS_STB="${IS_DB_DIR}/preterm_refseq_400.stb"
IS_GENES="${IS_DB_DIR}/preterm_refseq_400.genes.fna"

if [[ "$RUN_STEP6_INSTRAIN" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s6.done" ]]; then
    timestamp_start "inStrain"
    
    # 继承 Step 2 输出的双端 reads
    OUT1="${CLEAN_DIR}/${SAMPLE}.clean_1.fastq"
    OUT2="${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"
    
    # 定义输出路径 (改名为 preterm400)
    BAM_OUT="${QUANT_DIR}/${SAMPLE}_preterm400.sorted.bam"
    IS_OUT_DIR="${QUANT_DIR}/${SAMPLE}_inStrain.IS"

    # 2. Bowtie2 高质量比对与 BAM 流水线处理
    if [[ ! -f "$BAM_OUT" ]]; then
        echo "[inStrain] Running Bowtie2 mapping to preterm_refseq_400 database..."
        # 核心优化：--no-unal 不输出未比对上的 reads，大幅减小 BAM 文件体积
        bowtie2 -p "$THREADS" -x "$IS_BT2" -1 "$OUT1" -2 "$OUT2" --no-unal | \
        samtools view -bS - | \
        samtools sort -@ "$THREADS" -o "$BAM_OUT"
        
        echo "[inStrain] Indexing BAM file..."
        samtools index -@ "$THREADS" "$BAM_OUT"
    fi

    # 3. 运行 inStrain profile (火力全开模式)
    if [[ -f "$BAM_OUT" ]]; then
        echo "[inStrain] Running inStrain profile with full parameters..."
        
        # 挂载所有可用文件，开启 --database_mode 极速跳过未检出物种
        /data/wangjiaxuan/miniconda3/envs/meta1/bin/inStrain  \
        profile "$BAM_OUT" "$IS_FASTA" \
            -o "$IS_OUT_DIR" \
            -p "$THREADS" \
            -g "$IS_GENES" \
            -s "$IS_STB" \
            --database_mode
        
        # 4. 成功判定
        if [[ -d "${IS_OUT_DIR}/output" ]]; then
            touch "$LOG_DIR/${SAMPLE}.s6.done"
            timestamp_end
            
            # 5. inStrain 结果瘦身：删除无用的画图和臃肿的日志，保留 output 和 raw_data
            if [[ -d "${QUANT_DIR}/${SAMPLE}_inStrain.IS" ]]; then
                echo "Cleaning up inStrain bloated logs and figures..."
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
# 前提：依赖 Step 2 的标准双端 reads ($OUT1 和 $OUT2)

if [[ "$RUN_STEP7_SGVF" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s7.done" ]]; then
    timestamp_start "SGVFinder2"
    
    # 1. 绝对路径与环境定义
    ICRA_BIN="/data/wangjiaxuan/miniconda3/envs/sgvfinder/bin/icra"
    SVF_BIN="/data/wangjiaxuan/miniconda3/envs/sgvfinder/bin/svfinder"

    # 【核心修改 1：完美的统一前缀】指向你刚刚建好的新目录和新前缀
    SGVF_DB_PREFIX="/data/wangjiaxuan/refer/preterm_refseq_400/sgvfinder/preterm_refseq_400_sgvfinder"

    SGVF_OUT_DIR="${QUANT_DIR}/${SAMPLE}_SGVF"
    mkdir -p "$SGVF_OUT_DIR"
    
    OUT1="${CLEAN_DIR}/${SAMPLE}.clean_1.fastq"
    OUT2="${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"

    # 用 find 命令雷达扫描：检查是否已经跑完生成了 .jsdel 文件
    JSDEL_FILE=$(find "$SGVF_OUT_DIR" -maxdepth 1 -name "*.jsdel" | head -n 1)

    # 2. 运行 ICRA：使用正确的长参数
    if [[ -z "$JSDEL_FILE" ]]; then
        echo "[SGVFinder2] Running ICRA (Iterative Coverage-based Read Assignment)..."
        
        # 【极其关键的防呆】：强制清理上次失败可能遗留的坏 bam，防止 icra 误判跳过比对！
        rm -f "$SGVF_OUT_DIR"/*.bam
        
        $ICRA_BIN \
            --fq1 "$OUT1" \
            --fq2 "$OUT2" \
            --db "$SGVF_DB_PREFIX" \
            --outfol "$SGVF_OUT_DIR" \
            --threads "$THREADS"
            
        # 跑完之后，再次扫描获取最新生成的 .jsdel 文件绝对路径
        JSDEL_FILE=$(find "$SGVF_OUT_DIR" -maxdepth 1 -name "*.jsdel" | head -n 1)
    fi

    # 3. 运行 svfinder get_sample_map
    if [[ -n "$JSDEL_FILE" ]]; then
        echo "[SGVFinder2] Generating sample coverage map from: $(basename "$JSDEL_FILE")"
        
        # 【核心修改 2】：--db_path 也必须使用带前缀的完整路径！
        # 因为底层的 Python 逻辑是拿着这个前缀去贴 '.dlen' 后缀的
        $SVF_BIN get_sample_map --delta_file "$JSDEL_FILE" --db_path "$SGVF_DB_PREFIX"
        
        # 显式捕获刚才那条命令的退出状态码
        if [[ $? -eq 0 ]]; then
            # 状态码为 0，说明成功了！发通行证！
            touch "$LOG_DIR/${SAMPLE}.s7.done"
            rm -f "$SGVF_OUT_DIR"/*.bam "$SGVF_OUT_DIR"/*.pmp "$SGVF_OUT_DIR"/*.jspi "$SGVF_OUT_DIR"/*.jsdel
            timestamp_end
        else
            # 状态码非 0，说明报错了
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
    # 动态检查：默认可以清理。但只要有一个开启的模块没跑完，就不清理。
    CAN_CLEAN="true"
    [[ "$RUN_STEP3_MPA" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s3.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP4_HUMANN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP5_ARGS" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP6_INSTRAIN" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s6.done" ]] && CAN_CLEAN="false"
    [[ "$RUN_STEP7_SGVF" == "true" && ! -f "$LOG_DIR/${SAMPLE}.s7.done" ]] && CAN_CLEAN="false" # <-- 新增这一行

    if [[ "$CAN_CLEAN" == "true" ]]; then
        echo "[$(get_time)] >>> ALL TARGET STEPS SUCCESSFUL: $SAMPLE. Starting cleanup."
        # 1. 恢复：删除占用空间巨大的 SAM 临时文件 (MetaPhlAn 产生的)
        [[ -f "${QUANT_DIR}/${SAMPLE}.sam.bz2" ]] && rm -f "${QUANT_DIR}/${SAMPLE}.sam.bz2"

        # 2. 删除拼接版的 clean.fastq 以释放空间
        if [[ -f "${CLEAN_DIR}/${SAMPLE}.clean.fastq" ]]; then
            rm -f "${CLEAN_DIR}/${SAMPLE}.clean.fastq"
        fi

        # 3. 压缩双端 clean data 留档 (pigz 压缩后会自动删除原文件)
        for fq in "${CLEAN_DIR}/${SAMPLE}.clean_1.fastq" "${CLEAN_DIR}/${SAMPLE}.clean_2.fastq"; do
            if [[ -f "$fq" ]]; then
                echo "Compressing $fq with pigz..."
                pigz -f -p 8 "$fq"
            fi
        done
        
        # 4. 生成大功告成标记
        touch "$LOG_DIR/${SAMPLE}.all_done"
    else
        echo "[$(get_time)] >>> SKIP CLEANUP: $SAMPLE has incomplete steps (or failed). Check logs."
    fi
fi

# ==========================================
# 无论是否 Cleanup，正常跑完后主动清空当前样本的内存盘
[[ -d "$SHM_DIR" ]] && rm -rf "$SHM_DIR"
# ==========================================

# 脚本最末尾的结束提示
echo "======= [$(get_time)] FINISHED: $SAMPLE ======="
