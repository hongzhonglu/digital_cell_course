# -*- coding: utf-8 -*-
"""
02_AlphaGEM构建流程演示.py
==========================
课程7&8：非模式物种 GEM 自动构建 —— AlphaGEM 流程讲解

AlphaGEM (hongzhonglu/AlphaGEM) is a deep-learning pipeline that uses protein
language models (e.g. ESM) for homology annotation and automatic GEM
reconstruction of non-model yeasts. Unlike gapseq (an alignment-driven
R/shell pipeline), AlphaGEM uses pretrained protein language model embeddings
to detect remote homologs — especially effective for sparsely annotated
non-model yeasts.

This script is a "pipeline walk-through demo": real commands are printed as
strings (dry-run), with commented-out subprocess invocations attached; the
deep-learning reconstruction is not actually executed in class.
"""

print("=" * 70)
print("Step 0: AlphaGEM installation and environment setup (Linux + GPU preferred)")
print("=" * 70)
install_cmds = [
    # 克隆课程组开源的 AlphaGEM 仓库
    "git clone https://github.com/hongzhonglu/AlphaGEM.git",
    "cd AlphaGEM",
    # 建立 conda 环境（依赖 python、pytorch、transformers、fair-esm 等）
    "conda env create -f environment.yml",
    "conda activate AlphaGEM",
]
for cmd in install_cmds:
    print("  $", cmd)

print("\n" + "=" * 70)
print("Step 1: input — protein fasta of the non-model yeast")
print("=" * 70)
# 与 gapseq 类似，输入是目标物种的蛋白序列；区别在于后续注释手段
print("  Example: data/protein_fasta/Saccharomycodes_ludwigii.faa")
print("  Note: each protein is encoded into an embedding vector by a language")
print("        model such as ESM, then searched against protein embeddings of a")
print("        reference species (e.g. S. cerevisiae), so that reactions/GPRs of")
print("        the reference GEM (e.g. Yeast8/yeast-GEM) are mapped onto the target.")

print("\n" + "=" * 70)
print("Step 2: protein language model embedding + homology annotation")
print("=" * 70)
cmd_annotate = "python scripts/run_homology_annotation.py --input data/protein_fasta/target.faa --out output/annotation"
print("  $", cmd_annotate)
print("  Key idea: embedding similarity > sequence identity,")
print("  capturing remote homologous enzymes that blast misses, improving")
print("  annotation coverage for non-model yeasts.")
print("  Output: mapping table (csv) between target genes and reference genes/reactions")

print("\n" + "=" * 70)
print("Step 3: assemble the draft model from the mapping results")
print("=" * 70)
cmd_build_draft = "python scripts/build_draft_model.py --mapping output/annotation/mapping.csv --template template/yeast8.xml --out output/draft_model.xml"
print("  $", cmd_build_draft)
print("  Method: using the reference yeast GEM as template, reactions whose genes")
print("  have homologs in the target species are kept, and the GPR rules are")
print("  rewritten with the target species' genes.")

print("\n" + "=" * 70)
print("Step 4: gap-filling and growth validation")
print("=" * 70)
cmd_gapfill = "python scripts/gapfill_model.py --model output/draft_model.xml --medium config/minimal_medium.csv --out output/final_model.xml"
print("  $", cmd_gapfill)
print("  If the draft cannot grow on minimal medium, missing reactions are added")
print("  from a reaction database at minimal cost (gap-filling) until the")
print("  minimal growth-rate requirement is met.")

print("\n" + "=" * 70)
print("Step 5: downstream analysis of the final model")
print("=" * 70)
print("  The final model (final_model.xml) can be loaded directly with COBRApy:")
print("      from cobra.io import read_sbml_model")
print("      model = read_sbml_model('output/final_model.xml')")
print("  Then FBA, substrate-utilization screening and essential-gene prediction")
print("  can be performed (see scripts 03/04 of this module).")

print("\n[Comparison summary]")
print("  gapseq  : blast alignment -> R/shell pipeline -> common for prokaryotic non-model species")
print("  AlphaGEM: protein language model embeddings -> deep-learning pipeline,")
print("            strong at remote homology; suited to sparsely annotated non-model yeasts")

# 若在服务器上真实执行（课堂演示中注释掉）：
# import subprocess
# subprocess.run(cmd_annotate, shell=True, check=True)
