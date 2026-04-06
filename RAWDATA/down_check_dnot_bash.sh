while read -r id; do
    # 定义 R1 和 R2 的路径
    R1="./${id}_1.fastq.gz"
    R2="./${id}_2.fastq.gz"
    
    # 检查文件是否存在
    if [[ ! -f "$R1" ]] || [[ ! -f "$R2" ]]; then
        echo "--- Sample $id is incomplete ---"
        [[ ! -f "$R1" ]] && echo "  [MISSING] $R1"
        [[ ! -f "$R2" ]] && echo "  [MISSING] $R2"
    fi
done < sra.list


# 如果集群安装了 BBMap/BBTools
reformat.sh in=SRR5405966.fastq.gz \
            out1=SRR5405966_1.fastq.gz \
            out2=SRR5405966_2.fastq.gz