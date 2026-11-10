# -*- coding: utf-8 -*-
"""
01_gapseq构建流程演示.py
========================
课程7&8：非模式物种基因组尺度代谢模型（GEM）自动构建 —— gapseq 流程讲解

gapseq (Zimmermann et al., 2021, Genome Biology) is an R/shell based pipeline
for automatic GEM reconstruction and gap-filling of non-model organisms.
It must run on Linux (depends on R, sequence alignment tools and an LP solver),
so it is NOT executed live on classroom Windows/Mac laptops. This script is a
"pipeline walk-through demo":
  - Chinese comments + English print statements explain each step of gapseq;
  - Real command-line commands are shown as strings (dry-run printing, plus
    commented-out subprocess invocations), so students can copy them verbatim
    onto a Linux server.

Reference: course repo codes_part1/5_gapseq_using.py and gapseq docs
https://github.com/jotech/gapseq
"""

import os  # 仅用于打印示例路径，不执行真实命令

# =====================================================================
# Step 0: gapseq installation (run once on a Linux server)
# =====================================================================
print("=" * 70)
print("Step 0: gapseq installation (Linux server, run once)")
print("=" * 70)
# gapseq 以源码方式安装，依赖 R(>=3.6) 及若干 Bioconductor/CRAN 包
install_cmds = [
    "git clone https://github.com/jotech/gapseq ~/gapseq",   # 克隆源码仓库
    "cd ~/gapseq",
    # 安装 R 依赖包（gapseq 自带脚本，含 data.table、sybil、getopt 等）
    "Rscript ./src/install_dependencies.R",
    # 可选：安装支持的求解器（cplex / glpk；glpk 免费可用）
    "sudo apt-get install glpk-utils libglpk-dev",
]
for cmd in install_cmds:
    print("  $", cmd)

# 若在服务器上，可用 subprocess 真实调用（课堂演示中注释掉）：
# import subprocess
# subprocess.run(["git", "clone", "https://github.com/jotech/gapseq"], check=True)

# =====================================================================
# Step 0b (quick route): gapseq doall —— 课程仓库示例的一步式命令
# =====================================================================
# 课程仓库 codes_part1/5_gapseq_using.py 使用的就是 doall 一步式命令：
#     ./gapseq doall toy/myb71.fna
# 它把下面第 2~6 步（find / find-transport / draft / medium / fill）串成一条流水线，
# 适合快速体验；若要控制每步参数，仍建议按下文逐步执行。
print("\n" + "=" * 70)
print("Step 0b: quick route — 'gapseq doall' (as in course repo codes_part1/5_gapseq_using.py)")
print("=" * 70)
# doall 一步完成 find → draft → fill 全流程
cmd_doall = "~/gapseq/gapseq doall mystrain.fna"
print("  $", cmd_doall)
# 也可以指定自定义成分明确培养基做缺口填充（课程示例：MM_glu.csv）
cmd_doall_medium = "~/gapseq/gapseq doall mystrain.fna dat/media/MM_glu.csv"
print("  $", cmd_doall_medium, "   # gap-fill against a defined medium")
print("  Note: 'doall' chains Steps 2-6 below into one pipeline;")
print("        the step-by-step route gives finer control over each parameter.")

# =====================================================================
# Step 1: input — genome/protein fasta of the non-model organism
# =====================================================================
print("\n" + "=" * 70)
print("Step 1: prepare input fasta (genomic nucleotide or protein sequence)")
print("=" * 70)
# gapseq 的输入可以是一个属/种的全基因组核苷酸 fasta（.fna）
# 或注释好的蛋白 fasta（.faa）。非模式物种常从 NCBI 下载：
print("  $ wget 'https://ftp.ncbi.nlm.nih.gov/genomes/.../GCF_xxxx.fna.gz'")
print("  $ gunzip GCF_xxxx.fna.gz && mv GCF_xxxx.fna mystrain.fna")
print("  Note: gapseq internally uses tblastn / blastp to align sequences")
print("        against its own reaction-sequence database (MetaCyc/KEGG/ModelSEED).")

# =====================================================================
# Step 2: gapseq find —— 基于同源性的反应与通路注释
# =====================================================================
print("\n" + "=" * 70)
print("Step 2: gapseq find — pathway prediction and reaction annotation")
print("=" * 70)
# find 子命令对每条候选反应给出证据打分（bitscore 阈值 + 通路完整性）
cmd_find = "~/gapseq/gapseq find -p all -t mystrain.fna > mystrain_find.log 2>&1"
print("  $", cmd_find)
print("  Output: mystrain-Reactions.tbl (reaction table) and mystrain-Pathways.tbl (pathway table)")
print("  -p all   predict all pathways (can also restrict to core metabolism)")
print("  -t       input is nucleotide sequence -> tblastn; omit -t for protein fasta")

# =====================================================================
# Step 3: gapseq find-transport —— 转运反应预测
# =====================================================================
print("\n" + "=" * 70)
print("Step 3: gapseq find-transport — transporter prediction")
print("=" * 70)
cmd_transport = "~/gapseq/gapseq find-transport mystrain.fna"
print("  $", cmd_transport)
print("  Output: mystrain-Transporter.tbl (TCDB-based transporter annotation)")

# =====================================================================
# Step 4: gapseq draft —— 生成草稿模型（draft GEM）
# =====================================================================
print("\n" + "=" * 70)
print("Step 4: gapseq draft — assemble the draft genome-scale model")
print("=" * 70)
cmd_draft = ("~/gapseq/gapseq draft -r mystrain-Reactions.tbl "
             "-t mystrain-Transporter.tbl -p mystrain-Pathways.tbl "
             "-c mystrain.fna -u 200 -l 100")
print("  $", cmd_draft)
print("  -u 200  high-confidence bitscore cutoff for universal reactions (lower -> more reactions)")
print("  -l 100  low-confidence bitscore cutoff")
print("  Output: mystrain-draft.RDS / mystrain.xml (SBML draft model)")
print("  The draft cannot grow on the medium yet — gaps remain!")

# =====================================================================
# Step 5: gapseq medium —— 预测培养基成分
# =====================================================================
print("\n" + "=" * 70)
print("Step 5: gapseq medium — predict the required growth medium")
print("=" * 70)
cmd_medium = "~/gapseq/gapseq medium -m mystrain-draft.RDS -p mystrain-Pathways.tbl"
print("  $", cmd_medium)
print("  Output: mystrain-medium.csv (default medium derived from biomass precursor demands)")

# =====================================================================
# Step 6: gapseq fill —— 缺口填充，得到可生长的最终模型
# =====================================================================
print("\n" + "=" * 70)
print("Step 6: gapseq fill — gap-filling against that medium")
print("=" * 70)
cmd_fill = ("~/gapseq/gapseq fill -m mystrain-draft.RDS "
            "-n mystrain-medium.csv -c mystrain-Reactions.tbl "
            "-g 0.05 -b 100")
print("  $", cmd_fill)
print("  -g 0.05  required minimal growth rate; unmet -> add fewest reactions to close gaps")
print("  -b 100   evidence bitscore cutoff for added reactions")
print("  Output: mystrain.RDS and mystrain.xml — the final growth-capable GEM!")

# =====================================================================
# Step 7: 产物模型的质量评估与下游分析
# =====================================================================
print("\n" + "=" * 70)
print("Step 7: quality assessment and downstream analysis of the final model")
print("=" * 70)
print("  1) Quality: run MEMOTE for a model quality report (stoichiometric consistency etc.)")
print("     $ memote run mystrain.xml --filename memote_report.html")
print("  2) Downstream: load mystrain.xml into COBRApy for FBA, substrate-utilization")
print("     screening, essential-gene prediction (see scripts 03/04 of this module).")
print("  3) Quick in-R check:")
print("     Rscript -e 'library(sybil); m <- readRDS(\"mystrain.RDS\"); print(checkGrowth(m))'")

print("\n[Summary] fasta -> find -> find-transport -> draft -> medium -> fill")
print("turns a non-model organism's genome into a simulatable metabolic model.")
