# High-Speed Metagenomics Pipeline (v3.8)

这是一个专为 **大规模宏基因组数据 (10,000+ 样本)** 设计的自动化分析流水线。针对 **机械硬盘 (HDD)** 的并发 I/O 瓶颈进行了深度优化，具备工业级的稳健性、断点续传能力和自动清理机制。

## 🌟 核心特性

- **I/O 敏感型优化**：自动处理 HUMAnN 数据库组装过程中的磁头寻道冲突，平衡并发数与处理速度。
- **全流程断点续传**：每个步骤独立生成 `.done` 标记，支持随时中断并安全重启。
- **MetaPhlAn 4.2 适配**：修正了 Python 3.13 环境下的管道解析 Bug，强制数据落地以确保稳定性。
- **极致清洁机制**：引入 `trap` 信号捕获，任务崩溃或中止时自动物理抹除带随机后缀的 `humann_temp` 目录。
- **智能收尾清理**：只有在物种、功能、抗性基因三大核心分析全部成功后，才会触发 SAM 清理和 Fastq 压缩。

## 🛠️ 快速启动

### 第一步：构建元数据配置表
在开始批处理之前，需要扫描原始数据目录并生成任务清单。

```bash
python step1_build_check_config.py -i /path/to/your/rawdata

-i: 指定原始 Fastq 数据所在的文件夹。

其他可选参数请通过 python step1_build_check_config.py --help 查看。
```

### 第二步：启动并行流水线
确认配置表无误后，启动主调度脚本。该脚本默认开启 6 个样本并发。

```Bash
bash step2_run_pipeline.sh
```

### 📂 目录结构说明

```Plaintext
.
├── 1.Clean_Data/      # 质控与去污染后的结果 (Pigz 压缩保存)
├── 2.Quantity/        # 物种图谱 (MetaPhlAn) 与功能定量 (HUMAnN)
├── 3.Annotation/ARG/  # 耐药基因分析结果 (ARGs-OAP)
├── logs/              # 运行日志、Done 标记及 Runtime Joblog
├── Rawdata/           # 原始输入数据目录
├── step1_build_check_config.py # 环境预检与表构建
└── step2_run_pipeline.sh       # 多样本并发调度器
└── single_sample_process.sh    # 单样本核心逻辑脚本
```

### 🔧 模块流程

`centrifuger=1.1.0`: 极速去除人类宿主序列污染（基于 CHM13 索引）。

`fastp=1.1.0`: 自动质控、接头剪切并进行 Reads 合并。

`metaphlan=4.2.4`: 基于 SGB 数据库的物种组成定量，修正管道解析错误。

`HUMAnN=4.0.0a2`: 功能基因家族与通路定量（已优化磁头负载，移除了 --bypass-prescreen 参数以平衡组装开销）。

`ARGs-OAP=3.2.4`: 抗性基因定量，通过独立 work_space 规避并发文件名冲突。

## 📧 联系方式 (Contact)

如果您在使用本流水线过程中遇到任何问题，或者有合作研究的意向，欢迎通过以下方式联系：

- **作者**: [王家轩]
- **Email**: [poomourse@126.com]
- **GitHub**: [https://github.com/wangjiaxuan666](https://github.com/wangjiaxuan666)
- **个人网页**: [abego.cn](https://abego.cn/)

对于程序 Bug 或功能建议，建议直接在本项目仓库中提交 **Issue**。

---

## 更新计划

看项目需求。。。
