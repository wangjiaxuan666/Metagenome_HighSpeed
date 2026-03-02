#!/usr/bin/python3
import os
import sys
import argparse
import glob
from pathlib import Path

def get_args():
    parser = argparse.ArgumentParser(description="宏基因组大规模样本并行管理工具")
    subparsers = parser.add_subparsers(dest="command", help="执行模式: build (生成配置) 或 check (检查状态)")

    # Build 模式参数
    build_parser = subparsers.add_parser("build", help="扫描目录生成样本清单")
    build_parser.add_argument("-i", "--input", required=True, help="测序原始数据目录 (Raw Data)")
    build_parser.add_argument("-o", "--output", default="fq.meta", help="生成的 Meta 文件名 (默认: fq.meta)")
    build_parser.add_argument("--r1", default=".raw.R1.fq.gz", help="R1 文件的后缀 (默认: _R1.fq.gz)")
    build_parser.add_argument("--r2", default=".raw.R2.fq.gz", help="R2 文件的后缀 (默认: _R2.fq.gz)")

    # Check 模式参数
    check_parser = subparsers.add_parser("check", help="对比日志，提取未完成样本")
    check_parser.add_argument("-m", "--meta", required=True, help="原始的 Meta 配置文件")
    check_parser.add_argument("-l", "--log", default="logs", help="日志存放目录 (默认: logs)")
    check_parser.add_argument("-r", "--retry", default="retry.fq.meta", help="生成的补跑清单文件名")
    check_parser.add_argument("--flag", default="JOB FINISHED:", help="日志中的成功标志字符串")

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)
    return parser.parse_args()

def build_meta(args):
    """扫描目录并生成 3 列格式的 TSV"""
    input_dir = Path(args.input)
    r1_suffix = args.r1
    r2_suffix = args.r2
    
    # 使用 pathlib 提升路径处理兼容性
    r1_files = sorted(list(input_dir.glob(f"*{r1_suffix}")))
    
    if not r1_files:
        print(f"⚠️ 警告: 在 {input_dir} 中未找到匹配 {r1_suffix} 的文件。")
        return

    count = 0
    with open(args.output, 'w') as f:
        for r1_path in r1_files:
            sample_id = r1_path.name.replace(r1_suffix, "")
            # 自动推导 R2 路径
            r2_path = Path(str(r1_path).replace(r1_suffix, r2_suffix))
            
            if r2_path.exists():
                f.write(f"{sample_id}\t{r1_path.absolute()}\t{r2_path.absolute()}\n")
                count += 1
            else:
                print(f"❓ 样本 {sample_id} 缺少 R2 文件: {r2_path.name}")
                
    print(f"📋 清单已生成: {args.output} (共 {count} 个样本)")

def check_status(args):
    """高效对比日志"""
    if not os.path.exists(args.meta):
        print(f"❌ 错误: 找不到配置文件 {args.meta}")
        return

    # 1. 加载原始 Meta 数据
    all_samples = {}
    with open(args.meta, 'r') as f:
        for line in f:
            parts = line.strip().split('\t')
            if len(parts) >= 1:
                all_samples[parts[0]] = line.strip()

    # 2. 扫描日志 (针对上万样本，建议日志命名包含样本 ID 以便快速定位)
    success_ids = set()
    log_dir = Path(args.log)
    if not log_dir.exists():
        print(f"❌ 错误: 日志目录 {args.log} 不存在")
        return

    print(f"🔍 正在扫描 {args.log} 中的日志文件...")
    # 仅读取最后几行提高效率，或全量扫描
    for log_file in log_dir.glob("*.out"):
        try:
            with open(log_file, 'r', errors='ignore') as f:
                # 倒序读取最后 10 行通常就能找到成功标志，极大节省大日志读取时间
                lines = f.readlines()[-10:] 
                for line in lines:
                    if args.flag in line:
                        s_id = line.split(args.flag)[-1].strip()
                        success_ids.add(s_id)
        except Exception:
            continue

    # 3. 结果输出
    missing_ids = set(all_samples.keys()) - success_ids
    
    if not missing_ids:
        print("✅ 所有样本已处理完成！")
    else:
        with open(args.retry, 'w') as f:
            for m_id in sorted(list(missing_ids)):
                f.write(all_samples[m_id] + "\n")
        
        print(f"📊 统计: 总计 {len(all_samples)} | 成功 {len(success_ids)} | 缺失 {len(missing_ids)}")
        print(f"🛠️ 补跑清单: {args.retry}")
        print(f"💡 SLURM 提示: --array=1-{len(missing_ids)}")

if __name__ == "__main__":
    args = get_args()
    if args.command == "build":
        build_meta(args)
    elif args.command == "check":
        check_status(args)