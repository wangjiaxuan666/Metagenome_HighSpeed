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
RUN_STEP2_FASTP=${RUN_STEP2_FASTP:-true} # 核心替换：fastp
RUN_STEP3_MPA=${RUN_STEP3_MPA:-true}
RUN_STEP35_STRAIN=${RUN_STEP35_STRAIN:-true}
RUN_STEP4_HUMANN=${RUN_STEP4_HUMANN:-true}
RUN_STEP5_ARGS=${RUN_STEP5_ARGS:-true}
RUN_CLEANUP=${RUN_CLEANUP:-true}

# ---------- 1. 路径与环境配置 ----------
CONDA_BASE=$(conda info --base)
source "$CONDA_BASE/etc/profile.d/conda.sh"
conda activate meta

CENT_BIN="/data/wangjiaxuan/miniconda3/envs/meta/bin/centrifuger"
CENT_DB="/data/wangjiaxuan/refer/centrifuger_db/CHM13/chm13_index"
DB_PKL="/data/wangjiaxuan/refer/metaphlan_db/mpa_vOct22_CHOCOPhlAnSGB_202403.pkl"

CLEAN_DIR="1.Clean_Data"
QUANT_DIR="2.Quantity"
ARG_DIR="3.Annotation/ARG"
LOG_DIR="logs"
mkdir -p $CLEAN_DIR $QUANT_DIR $ARG_DIR $LOG_DIR

LOG_FILE="$LOG_DIR/${SAMPLE}.process.log"
exec > >(tee -a "$LOG_FILE") 2>&1

# ---------- 2. 工具函数 ----------

function timestamp_start() {
    MODULE_NAME=$1
    CUR_START=$(date +%s)
    echo "[$(date)] >>> START: $MODULE_NAME"
}

function timestamp_end() {
    CUR_END=$(date +%s)
    DIFF=$((CUR_END - CUR_START))
    echo "[$(date)] <<< END: $MODULE_NAME | Time: ${DIFF}s"
    echo "[TIME_AUDIT] Sample: $SAMPLE | Module: $MODULE_NAME | Time: ${DIFF}s" >> $LOG_DIR/time_summary.audit
}

function check_disk_safety() {
    FREE_MB=$(df -m /data | tail -1 | awk '{print $4}')
    [[ "$FREE_MB" -lt 153600 ]] && { echo "Disk space low, waiting 10m..."; sleep 600; }
}

# ---------- 3. 执行流程 ----------

echo "======= [$(date)] Processing: $SAMPLE ======="

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

# --- Step 2: fastp (质控 - 仅保留 JSON) ---
if [[ "$RUN_STEP2_FASTP" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s2.done" ]]; then
    timestamp_start "fastp"
    IN1="${CLEAN_DIR}/${SAMPLE}_nonhuman_1.fq.gz"
    IN2="${CLEAN_DIR}/${SAMPLE}_nonhuman_2.fq.gz"
    OUT_FQ="${CLEAN_DIR}/${SAMPLE}.clean.fastq"

    # 增加 -V (suppress html report) 选项
    fastp -i "$IN1" -I "$IN2" \
          -o /dev/null -O /dev/./null \
          --merged_out "$OUT_FQ" \
          --merge \
          -q 20 -u 30 -l 50 \
          -w "$THREADS" \
          -j "${LOG_DIR}/${SAMPLE}.fastp.json" \
          -h /dev/null \
          --thread "$THREADS"

    if [[ -f "$OUT_FQ" ]]; then
        rm -f "$IN1" "$IN2"
        touch "$LOG_DIR/${SAMPLE}.s2.done"
        timestamp_end
    fi
fi

# 后续流程的输入文件
CLEAN_FQ="${CLEAN_DIR}/${SAMPLE}.clean.fastq"

# --- Step 3: MetaPhlAn 4 (带失败自动清理) ---
if [[ "$RUN_STEP3_MPA" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s3.done" ]]; then
    timestamp_start "MetaPhlAn"
    # 强制落地 bt2out 解决管道 Bug
    metaphlan "$CLEAN_FQ" -o "${QUANT_DIR}/${SAMPLE}_profile.txt" \
              --input_type fastq --nproc "$THREADS" --offline --sample_id "${SAMPLE}" \
              -x mpa_vOct22_CHOCOPhlAnSGB_202403 -s "${QUANT_DIR}/${SAMPLE}.sam.bz2" \
              --mapout "${QUANT_DIR}/${SAMPLE}.bt2out" \
              -t rel_ab_w_read_stats

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
    if [[ -f "${QUANT_DIR}/${SAMPLE}.sam.bz2" ]]; then
        sample2markers.py -i "${QUANT_DIR}/${SAMPLE}.sam.bz2" -o "${QUANT_DIR}/db_markers/" \
                          -n "$THREADS" -f bz2 -d "$DB_PKL" -b 80 \
                          --min_reads_aligning 8
        if [[ -f "${QUANT_DIR}/db_markers/${SAMPLE}.markers.pkl" ]]; then
            rm -f "${QUANT_DIR}/${SAMPLE}.sam.bz2"
            touch "$LOG_DIR/${SAMPLE}.s3_5.done"
            timestamp_end
        fi
    fi
fi

# --- Step 4: HUMAnN (针对随机后缀优化的清理版) ---
if [[ "$RUN_STEP4_HUMANN" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s4.done" ]]; then
    timestamp_start "HUMAnN"
    
    # 1. 设置清理陷阱 (Trap)
    # 只要任务异常退出，立即删除该样本名下所有的 humann_temp 文件夹
    # 使用 ${SAMPLE}_humann_temp* 匹配随机后缀
    trap 'echo "[CLEANUP] Interrupted! Wiping temp files..."; rm -rf "${QUANT_DIR}/${SAMPLE}_humann_temp"* ; exit 1' SIGINT SIGTERM ERR

    # 2. 执行 HUMAnN
    # --o-log: 记录日志，方便以后排查为什么跑得慢
    humann --threads "$THREADS" \
           --input "$CLEAN_FQ" \
           --output "$QUANT_DIR" \
           --output-basename "$SAMPLE" \
           --taxonomic-profile "${QUANT_DIR}/${SAMPLE}_profile.txt" \
           --remove-temp-output

    # 3. 结果验证与强制清理
    # 注意：检查标准输出文件名是否为 ${SAMPLE}_genefamilies.tsv
    if [[ -f "${QUANT_DIR}/${SAMPLE}_2_genefamilies.tsv" ]]; then
        trap - SIGINT SIGTERM ERR  # 成功后解除陷阱
        touch "$LOG_DIR/${SAMPLE}.s4.done"
        timestamp_end
    else
        echo "[ERROR] HUMAnN failed for $SAMPLE. Cleaning up..."
        rm -rf "${QUANT_DIR}/${SAMPLE}_humann_temp"*
        exit 1
    fi
fi

# --- Step 5: ARGs-OAP (带强制清理逻辑) ---
if [[ "$RUN_STEP5_ARGS" == "true" ]] && [[ ! -f "$LOG_DIR/${SAMPLE}.s5.done" ]]; then
    timestamp_start "ARGs_OAP"
    W_DIR="${ARG_DIR}/${SAMPLE}" ; mkdir -p "${W_DIR}/work_space"
    ln -sf $(readlink -f "$CLEAN_FQ") "${W_DIR}/work_space/sample.fastq"
    
    PUSH_DIR=$(pwd)
    cd "${W_DIR}/work_space" || exit

    args_oap stage_one -i ./ -o ./result -f fastq -t "$THREADS"
    args_oap stage_two -i ./result -t "$THREADS" -o ./result
    
    if [[ -d "./result" ]]; then
        mv ./result/normalized_cell.*.txt ../
        cd "$PUSH_DIR" || exit
        rm -rf "${W_DIR}/work_space" # 成功清理
        [[ -f "${PUSH_DIR}/fifo_map.mapout.txt" ]] && rm -f "${PUSH_DIR}/fifo_map.mapout.txt"
        touch "$LOG_DIR/${SAMPLE}.s5.done" ; timestamp_end
    else
        cd "$PUSH_DIR" || exit
        echo "[ERROR] ARGs-OAP failed. Wiping work_space..."
        rm -rf "${W_DIR}/work_space" # 失败也强制清理，不留垃圾
        [[ -f "${PUSH_DIR}/fifo_map.mapout.txt" ]] && rm -f "${PUSH_DIR}/fifo_map.mapout.txt"
        exit 1
    fi
fi

# --- Final Cleanup (严谨全通过版) ---
if [[ "$RUN_CLEANUP" == "true" ]]; then
    # 必须满足：S3(物种)、S4(功能)、S5(抗性基因) 全部有 .done 标记
    if [[ -f "$LOG_DIR/${SAMPLE}.s3.done" ]] && \
       [[ -f "$LOG_DIR/${SAMPLE}.s4.done" ]] && \
       [[ -f "$LOG_DIR/${SAMPLE}.s5.done" ]]; then
        
        echo "[$(date)] >>> ALL STEPS SUCCESSFUL: $SAMPLE. Starting cleanup."

        # 1. 删除占用空间巨大的 SAM 临时文件 (MetaPhlAn 产生的)
        [[ -f "${QUANT_DIR}/${SAMPLE}.sam.bz2" ]] && rm -f "${QUANT_DIR}/${SAMPLE}.sam.bz2"
        
        # 2. 压缩 Fastq 结果 (增加判断，防止重复压缩)
        # 只有在 .fastq 存在且 .fastq.gz 不存在时才执行压缩
        if [[ -f "$CLEAN_FQ" ]] && [[ ! -f "${CLEAN_FQ}.gz" ]]; then
            echo "Compressing $CLEAN_FQ with pigz..."
            pigz -p 8 "$CLEAN_FQ"
        fi
        
        # 3. 生成一个最终的“大功告成”标记（可选，方便统计）
        touch "$LOG_DIR/${SAMPLE}.all_done"
    else
        # 只要有一个没过，就保持原样，方便重跑排查
        echo "[$(date)] >>> SKIP CLEANUP: $SAMPLE has incomplete steps. Check logs."
    fi
fi

echo "======= [$(date)] FINISHED: $SAMPLE ======="