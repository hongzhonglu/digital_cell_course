# -*- coding: utf-8 -*-
# =====================================================================
# 03_单基因敲除与必需基因预测.py
# ---------------------------------------------------------------------
# 模块二《GEMs构建与表型预测》（对应课程4）第 3 个上机脚本。
#
# 目的：使用 cobra.flux_analysis.single_gene_deletion 对 e_coli_core 做
#       全模型单基因敲除模拟：
#       1) 统计必需基因 / 非必需基因数量；
#       2) 列出生长受损最严重的基因；
#       3) 结果导出 CSV 并画“敲除后生长速率分布”直方图。
#
# 必需基因判定标准：敲除后最大生长速率 < 野生型的 1%（视为致死）。
#
# 用法：python 03_单基因敲除与必需基因预测.py
# 依赖：cobra、pandas、matplotlib
# =====================================================================
"""Genome-wide single-gene deletion simulation on e_coli_core to predict essential genes."""

import os

import matplotlib
matplotlib.use("Agg")                       # 无显示环境保存图片
import matplotlib.pyplot as plt
import pandas as pd

from cobra.flux_analysis import single_gene_deletion
from model_loader import load_ecoli_core

# 输出目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False

model = load_ecoli_core()
print("Model:", model.id, f"({len(model.genes)} genes)")

# ----------------------------------------------------------------------
# 1. 先求野生型（未敲除）最大生长速率作为参照
# ----------------------------------------------------------------------
wt_growth = model.optimize().objective_value       # 野生型生长速率
print(f"\nWild-type max growth rate: {wt_growth:.4f} 1/h")

# ----------------------------------------------------------------------
# 2. 全模型单基因敲除模拟
# ----------------------------------------------------------------------
# single_gene_deletion 会依次“敲除”每个基因（关闭其 GPR 控制的反应），
# 返回 DataFrame：ids(基因集合) + growth(敲除后生长速率) + status
print("\nRunning single-gene deletion simulation "
      f"({len(model.genes)} genes in total)...")
ko_result = single_gene_deletion(model)            # 运行敲除分析
print("Simulation finished, result rows:", len(ko_result))
print("\nFirst 5 rows of the result table:")
print(ko_result.head().to_string(index=False))

# ----------------------------------------------------------------------
# 3. 整理结果：提取基因 id，计算相对生长速率，判定必需性
# ----------------------------------------------------------------------
# ko_result["ids"] 是集合类型（单基因敲除时每个集合只有 1 个基因）
ko_result["gene"] = ko_result["ids"].apply(lambda s: list(s)[0])
# 不可行解时 growth 为 NaN，fillna(0) 记为不能生长
ko_result["growth"] = ko_result["growth"].fillna(0.0)
# 相对生长速率 = 敲除后 / 野生型
ko_result["relative_growth"] = ko_result["growth"] / wt_growth

ESSENTIAL_CUTOFF = 0.01                            # 必需基因阈值：相对生长 < 1%
ko_result["is_essential"] = ko_result["relative_growth"] < ESSENTIAL_CUTOFF

essential = ko_result[ko_result["is_essential"]]     # 必需基因子表
nonessential = ko_result[~ko_result["is_essential"]] # 非必需基因子表

print("\n===== Essential gene statistics =====")
print(f"Total genes: {len(ko_result)}")
print(f"Essential genes (knock-out lethal): {len(essential)}")
print(f"Non-essential genes: {len(nonessential)}")
print(f"Essential fraction: {len(essential)/len(ko_result)*100:.1f}%")

# ----------------------------------------------------------------------
# 4. 展示前 20 个必需基因（按敲除后生长速率升序）
# ----------------------------------------------------------------------
print("\nTop 20 essential genes (smaller growth = more essential):")
top_essential = essential.sort_values("growth").head(20)
for _, row in top_essential.iterrows():
    print(f"  {row['gene']:8s} growth after KO = {row['growth']:.4f} 1/h "
          f"(relative {row['relative_growth']*100:.2f}%)")

# 展示几个“部分必需”基因：生长显著下降但未致死（1% ~ 50%）
partial = ko_result[(ko_result["relative_growth"] >= ESSENTIAL_CUTOFF)
                    & (ko_result["relative_growth"] < 0.5)]
print(f"\nPartially essential genes (relative growth 1%~50%, {len(partial)} in total):")
for _, row in partial.sort_values("growth").iterrows():
    print(f"  {row['gene']:8s} relative growth after KO = {row['relative_growth']*100:.1f}%")

# ----------------------------------------------------------------------
# 5. 结果导出 CSV
# ----------------------------------------------------------------------
out = ko_result[["gene", "growth", "relative_growth", "is_essential", "status"]]
out = out.sort_values("growth")                    # 按敲除后生长速率排序
csv_path = os.path.join(OUTPUT_DIR, "单基因敲除结果.csv")
out.to_csv(csv_path, index=False, encoding="utf-8-sig")
print("\nFull knock-out results saved to:", csv_path)

# ----------------------------------------------------------------------
# 6. 画直方图：敲除后相对生长速率的分布
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.hist(ko_result["relative_growth"], bins=30, color="#2ca02c",
        edgecolor="black", alpha=0.8)
ax.axvline(ESSENTIAL_CUTOFF, color="red", linestyle="--",
           label=f"Essential cutoff = {ESSENTIAL_CUTOFF}")
ax.set_xlabel("Relative growth after knockout (mutant / WT)")
ax.set_ylabel("Number of genes")
ax.set_title("Single gene deletion: growth distribution (e_coli_core)")
ax.legend()
ax.grid(alpha=0.3)
fig_path = os.path.join(OUTPUT_DIR, "单基因敲除_生长分布直方图.png")
fig.savefig(fig_path, dpi=150, bbox_inches="tight")
print("Histogram saved to:", fig_path)

print("\n=== Script 03_单基因敲除与必需基因预测.py finished ===")
