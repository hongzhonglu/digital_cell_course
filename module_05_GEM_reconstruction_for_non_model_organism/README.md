# 模块05：非模式物种 GEMs 自动构建与分析

**对应课次：第 7 & 8 讲**

## 目标
1. 掌握 gapseq（R/shell，序列比对驱动）自动构建 GEM 的完整流程：fasta → find → find-transport → draft → medium → fill，以及课程仓库示例的一步式命令 `gapseq doall`。
2. 掌握 AlphaGEM（深度学习/蛋白质语言模型同源注释，hongzhonglu/AlphaGEM）构建非模式酵母 GEM 的流程。
3. 对自动构建产物模型做两类典型分析：底物利用筛选、必需基因预测。

## 依赖
```bash
pip install cobra pandas matplotlib requests --proxy http://10.86.13.73:5900
```
Python 3；gapseq/AlphaGEM 本身不安装（课堂演示流程，命令以打印形式给出）。

## 运行顺序与预期输出
| 脚本 | 内容 | 预期输出 |
|---|---|---|
| `01_gapseq构建流程演示.py` | 打印讲解 gapseq 各步骤真实命令（dry-run，可照抄到 Linux）；含课程仓库 `codes_part1/5_gapseq_using.py` 使用的 `gapseq doall` 一步式命令及其自定义培养基用法 | 终端流程说明（英文 print + 中文注释），无文件输出 |
| `02_AlphaGEM构建流程演示.py` | 打印讲解 AlphaGEM 流程与输入输出 | 终端流程说明，无文件输出 |
| `03_模型分析_底物利用.py` | 加载 iAF692（BiGG 下载，自动缓存；无网回退 e_coli_core），逐碳源 FBA | `substrate_utilization.csv`；终端列出可利用底物（iAF692 验证得到甲醇、丙酮酸、乙酸、甲胺类等 8 种可利用底物，符合其甲基营养菌特性） |
| `04_模型分析_必需基因预测.py` | single_gene_deletion 必需基因预测 | `essential_genes.csv`（iAF692 约 203/692 ≈ 29% 必需基因）、`all_gene_deletion.csv` |

## 代码风格说明（第二轮完善）
- 变量名/函数名为英文 snake_case；print 与 docstring 为英文；
- 解释性注释保留中文逐行/逐块；每个 `.py` 文件开头保留中文块注释（含课程定位）。

## 说明
- 03/04 脚本会优先读取同目录下的 `iAF692.json`（首次运行自动从 BiGG 下载）；下载失败时自动回退到 COBRApy 内置 e_coli_core 模型，保证任何环境均可复现算法流程。
- iAF692 是甲基营养菌 Methylobacterium extorquens AM1 的模型，存在烟酸/对氨基苯甲酸营养缺陷，03 脚本内含"必需补充成分"自动识别逻辑。
- gapseq 分步命令与 gapseq 官方文档（https://github.com/jotech/gapseq）一致；`doall` 快速路线对照课程仓库 `codes_part1/5_gapseq_using.py`。
