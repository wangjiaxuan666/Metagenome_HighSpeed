# 1. 清空旧清单

# 2. 遍历检查，自动删标记并生成新清单
for f in logs/*.s3.done; do
    sample=$(basename "$f" .s3.done)
    # 如果没有病毒丰度表
    if [ ! -f "2.Quantity/${sample}_vsc_profile.txt" ]; then
        echo "🚨 准备重置 $sample 以补跑 MetaPhlAn..."
        
        # 【关键】：必须同时删掉 s3 和全局完成标记！
        rm -f "logs/${sample}.s3.done" "logs/${sample}.all_done"
        
        # 从你的全局 fq.meta 里把这个样本的那一行抓出来，塞进新清单
        grep "^${sample}\b" fq.meta_all >> rerun_vsc.meta
    fi
done

echo "✅ 搞定！共找到 $(wc -l < rerun_vsc.meta) 个需要补跑的样本。"