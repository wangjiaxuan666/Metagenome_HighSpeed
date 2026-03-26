#!/bin/bash
# ==============================================================================
# 宏基因组流水线绝对安全收尾脚本 (Step 3: Secure Harvest)
# 特性：支持外部传参；采用 字节级比对 确保 cp 成功后才允许执行 rm
# ==============================================================================

# 1. 接收命令行外部参数
COHORT_NAME=$1

# 参数防呆检查：如果没有输入参数，直接报错并提示用法
if [[ -z "$COHORT_NAME" ]]; then
    echo "❌ 致命错误：未提供队列/项目名称！"
    echo "💡 正确用法：bash step3_secure_harvest.sh <队列名称>"
    echo "   例如操作：bash step3_secure_harvest.sh Project_Preterm_Batch1"
    exit 1
fi

# 建立该队列专属的 Clean Reads 下载转移目录
DOWNLOAD_DIR="${COHORT_NAME}_CleanReads_For_DL"
mkdir -p "$DOWNLOAD_DIR"

echo "=========================================================="
echo "🚀 开始执行 [ $COHORT_NAME ] 队列的安全转移任务..."
echo "=========================================================="

# ==============================================================================
# 🛡️ 核心护城河：安全复制并验证删除函数
# 逻辑：cp 前检查 -> 执行 cp -> 比对源文件与目标文件的字节大小 -> 完全一致才 rm
# ==============================================================================
safe_cp_then_rm() {
    local src="$1"
    local dest_dir="$2"
    local filename=$(basename "$src")
    local dest_file="${dest_dir}/${filename}"

    # 检查 1：cp 之前检查源文件是否存在
    if [[ ! -f "$src" ]]; then
        echo "    [跳过] 找不到源文件 (可能已转移): $filename"
        return 1
    fi

    echo "    -> 正在转移: $filename ..."
    # 执行复制
    cp "$src" "$dest_file"

    # 检查 2：cp 之后，核对文件是否生成，且【字节大小】是否 100% 一致
    local src_size=$(stat -c%s "$src")
    local dest_size=$(stat -c%s "$dest_file" 2>/dev/null || echo "0")

    if [[ -f "$dest_file" ]] && [[ "$src_size" -eq "$dest_size" ]]; then
        # 终极确认通过，安全删除源文件！
        rm -f "$src"
        echo "    [✅ 安全] 字节级校验通过 ($dest_size bytes)。已清理原文件。"
    else
        # 只要差一个字节，绝对不删原文件！
        echo "    [💥 警告] 复制异常或大小不符！(源:$src_size 目标:$dest_size)。保留原文件以防丢失！"
    fi
}
# ==============================================================================

# 2. 扫描所有成功跑到终点的样本 (以 .all_done 为金标准)
SUCCESS_SAMPLES=$(ls logs/*.all_done 2>/dev/null | awk -F'/' '{print $2}' | sed 's/\.all_done//')
COUNT=$(echo "$SUCCESS_SAMPLES" | wc -w)

if [[ $COUNT -eq 0 ]]; then
    echo "❌ 没找到任何带有 .all_done 标记的成功样本！"
    exit 1
fi

echo "✅ 共检测到 $COUNT 个完美通关的样本！开始执行字节级安全转移..."

# 3. 循环处理每个样本的 Clean Reads
for SAMPLE in $SUCCESS_SAMPLES; do
    echo "--------------------------------------------------"
    echo "📦 处理样本: $SAMPLE"
    
    FQ1="1.Clean_Data/${SAMPLE}.clean_1.fastq.gz"
    FQ2="1.Clean_Data/${SAMPLE}.clean_2.fastq.gz"
    
    # 调用安全函数转移双端数据
    safe_cp_then_rm "$FQ1" "$DOWNLOAD_DIR"
    safe_cp_then_rm "$FQ2" "$DOWNLOAD_DIR"
done

echo "=========================================================="
echo "🎉 阶段转移完毕！"
echo "📁 所有的 Clean Data 已安全存放在 [ $DOWNLOAD_DIR ] 文件夹中，随时可供下载。"
echo "📂 你的 2.Quantity, 3.Annotation 等所有分析结果原封不动地保留在原位！"
echo "=========================================================="

# ==============================================================================
# 🚨 [危险区] 可选：自动删除原始数据 (Raw Data)
# ==============================================================================
# 如果你确定想让脚本顺手把原始数据删掉，请取消下面这段代码的注释。
# 逻辑：它会去 fq.meta 里找到该样本的原始 R1/R2 路径，并在确认 Clean Data 安全落地后，删除原始数据。

echo "🗑️ 正在清理原始测序数据 (Raw Data)..."
for SAMPLE in $SUCCESS_SAMPLES; do
    # 确认 Clean Data 确实已经在下载文件夹里了，才敢去删原始数据
    if [[ -f "$DOWNLOAD_DIR/${SAMPLE}.clean_1.fastq.gz" ]]; then
        # 从 fq.meta 提取原始数据绝对路径
        RAW_R1=$(grep "^${SAMPLE}\b" fq.meta | awk '{print $2}')
        RAW_R2=$(grep "^${SAMPLE}\b" fq.meta | awk '{print $3}')
        
        [[ -f "$RAW_R1" ]] && rm -f "$RAW_R1" && echo "    [已删除] 原始数据 $RAW_R1"
        [[ -f "$RAW_R2" ]] && rm -f "$RAW_R2" && echo "    [已删除] 原始数据 $RAW_R2"
    fi
done

# echo "✅ 原始数据清理完毕！"