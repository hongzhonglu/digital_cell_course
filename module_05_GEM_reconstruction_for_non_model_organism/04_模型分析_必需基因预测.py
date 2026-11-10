# -*- coding: utf-8 -*-
"""
04_模型分析_必需基因预测.py
===========================
课程7&8：对"自动构建产物模型"做单基因敲除（single gene deletion）分析，
预测必需基因（essential genes）。

Essential-gene prediction for an automatically reconstructed GEM via
single-gene deletion analysis.

Model source is the same as in script 03: local/BiGG iAF692 first,
falling back to e_coli_core. A gene is called essential if its knockout
growth rate drops below 1% of wild type (or below an absolute threshold).
Results are saved to essential_genes.csv.
"""

import os
import pandas as pd
from cobra.io import load_json_model
from cobra.flux_analysis import single_gene_deletion  # COBRApy 单基因敲除工具
__file__='module_05_GEM_reconstruction_for_non_model_organism'
SCRIPT_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_JSON = os.path.join(SCRIPT_DIR, "iAF692.json")

# ---------------- 1. Load the model (same A/B strategy as script 03) ----------------
if os.path.exists(MODEL_JSON):
    model = load_json_model(MODEL_JSON)
    print(f"[Model] iAF692: {len(model.genes)} genes")
else:
    from cobra.io import load_model
    model = load_model("textbook")
    print(f"[Model] Fallback to e_coli_core: {len(model.genes)} genes")

try:
    model.solver = "glpk"   # 使用免费求解器
except Exception:
    pass

# ---------------- 2. Wild-type growth rate ----------------
wt_growth = model.optimize().objective_value   # 未敲除时的最大生长率
print(f"[Wild type] maximal growth rate = {wt_growth:.4f} 1/h")

# ---------------- 3. Single-gene deletion simulation ----------------
print("[Knockout] running single_gene_deletion over all genes ...")
# single_gene_deletion 返回 DataFrame：ids（基因集合）、growth、status
deletion_result = single_gene_deletion(model)

# ---------------- 4. Determine essential genes ----------------
ESSENTIAL_FRACTION = 0.01          # 敲除后生长率 < 野生型 1% 判定为必需
ABS_THRESHOLD = 1e-4               # 绝对阈值，防止 wt 极小时的数值问题

deletion_result["ids"] = deletion_result["ids"].apply(lambda s: list(s)[0])  # 集合取单基因
deletion_result["growth"] = deletion_result["growth"].fillna(0.0)            # 致死=NaN→0
deletion_result["essential"] = (deletion_result["growth"] < ESSENTIAL_FRACTION * wt_growth) | \
                               (deletion_result["growth"] < ABS_THRESHOLD)

essential_df = deletion_result[deletion_result["essential"]].copy()

# 附加基因名称信息（便于解读）
essential_df["gene_name"] = essential_df["ids"].map(
    {gene.id: gene.name for gene in model.genes})

# ---------------- 5. Save and print results ----------------
n_total = len(model.genes)
n_essential = len(essential_df)
print(f"\n[Result] essential genes: {n_essential} / {n_total} "
      f"({n_essential/n_total:.1%} of all genes)")
print(essential_df[["ids", "gene_name", "growth"]].head(20).to_string(index=False))
if n_essential > 20:
    print(f"... remaining {n_essential-20} genes are in the csv file")

out_csv = os.path.join(SCRIPT_DIR, "essential_genes.csv")
essential_df[["ids", "gene_name", "growth"]].rename(
    columns={"ids": "gene_id", "growth": "growth_after_knockout"}
).to_csv(out_csv, index=False, encoding="utf-8-sig")

# 同时保存全部敲除结果，便于画分布
deletion_result.rename(columns={"ids": "gene_id"})[
    ["gene_id", "growth", "essential"]
].to_csv(os.path.join(SCRIPT_DIR, "all_gene_deletion.csv"),
         index=False, encoding="utf-8-sig")

print(f"\nSaved: {out_csv}")
print("       all_gene_deletion.csv (growth after every gene knockout)")


