#!/bin/bash
# ==============================================================================
# 大规模并行调度器 - 生产环境版 (v3.1)
# ==============================================================================

META_FILE="fq.meta"
MAX_JOBS=5            # 64核服务器+机械硬盘，建议并发 4-6，防止 IO 阻塞
THREADS_PER_JOB=16    # 每个样本分配 10 线程 (6*10=60核)
JOBLOG="logs/parallel_runtime.joblog" # 核心：记录成败，支持断点续跑

# 确保 logs 目录存在
mkdir -p logs

echo "[$(date)] >>> 样本批次任务启动 <<<"
echo "并发任务数: $MAX_JOBS | 单任务线程: $THREADS_PER_JOB"
echo "------------------------------------------------"

# --- GNU Parallel 核心命令 ---
# --joblog: 极其重要！记录每个样本的退出状态
# --resume: 如果任务中断，下次运行会自动跳过已完成的样本
# --progress: 实时显示进度条和预计完成时间 (ETA)
# --timeout 18h: 保证时间不能超过8h

cat "$META_FILE" | /data/wangjiaxuan/miniconda3/envs/meta/bin/parallel \
    -j "$MAX_JOBS" \
    --timeout 18h \
    --colsep '\t' \
    --joblog "$JOBLOG" \
    --progress \
    "bash single_sample_process.sh {1} {2} {3} $THREADS_PER_JOB"

# 任务收尾统计
echo "------------------------------------------------"
# 直接数一下生成了多少个“大功告成”标记文件
SUCCESS_COUNT=$(ls logs/*.all_done 2>/dev/null | wc -l)
echo "[$(date)] 批次处理结束。"
echo "成功样本数: $SUCCESS_COUNT"
echo "详情请查阅: $JOBLOG"
