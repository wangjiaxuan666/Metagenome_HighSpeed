#!/bin/bash

MAX_RETRIES=10  # 最大重试次数
COUNT=0
SLEEP_TIME=10

while [ $COUNT -lt $MAX_RETRIES ]; do
    echo "开始/恢复下载任务 (第 $((COUNT+1)) 次尝试)..."
    
    /data/wangjiaxuan/miniforge3/bin/iseq -i sra.list -g -p 4 -t 4
    
    if [ $? -eq 0 ]; then
        echo "🎉 所有 SRA 数据下载完成！"
        exit 0
    else
        COUNT=$((COUNT+1))
        echo "⚠️ 第 $COUNT 次下载失败，还剩 $((MAX_RETRIES-COUNT)) 次机会，${SLEEP_TIME}秒后重试..."
        sleep $SLEEP_TIME
    fi
done

echo "❌ 已达到最大重试次数，部分数据可能下载失败，请人工检查 sra.list。"
exit 1