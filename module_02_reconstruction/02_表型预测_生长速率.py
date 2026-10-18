# -*- coding: utf-8 -*-
# =====================================================================
# 02_表型预测_生长速率.py
# ---------------------------------------------------------------------
# 模块二《GEMs构建与表型预测》（对应课程4）第 2 个上机脚本。
#
# 目的：以 e_coli_core 为例进行两类表型预测：
#       1) 不同碳源（葡萄糖/果糖/富马酸/乙酸/乳酸/琥珀酸/丙酮酸）下的最大生长速率；
#       2) 葡萄糖摄取速率扫描 —— 简化版“表型相平面(PhPP)”分析，
#          即鲁棒性分析(robustness analysis)：x 轴为葡萄糖摄取速率，
#          y 轴为最大生长速率，观察有氧/无氧条件下的差异。
#       结果存 CSV，并画 matplotlib 图存 PNG。
#
# 用法：python 02_表型预测_生长速率.py
# 依赖：cobra、pandas、matplotlib
# =====================================================================
"""Predict growth rates on different carbon sources and run a glucose-uptake
robustness (phenotype phase plane) scan with e_coli_core."""

import os

import matplotlib
matplotlib.use("Agg")                       # 无显示环境下用 Agg 后端保存图片
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd

from model_loader import load_ecoli_core

# 输出目录：与本脚本同级的 output/ 文件夹
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# 设置字体与负号显示（图中文字为英文，DejaVu Sans 即可正常显示）
plt.rcParams["font.sans-serif"] = ["DejaVu Sans"]
plt.rcParams["axes.unicode_minus"] = False  # 正常显示负号

model = load_ecoli_core()
print("Model:", model.id, f"({len(model.reactions)} reactions)")

# ----------------------------------------------------------------------
# 1. 不同碳源下的最大生长速率预测
# ----------------------------------------------------------------------
# e_coli_core 中可用的碳源交换反应（代谢物 id 可在模型中查到）
carbon_sources = {
    "Glucose":    "EX_glc__D_e",
    "Fructose":   "EX_fru_e",
    # 注意：e_coli_core 是核心模型，不含甘油交换反应（EX_glyc_e），
    # 若换用 iJO1366 等全基因组模型即可加入 "Glycerol": "EX_glyc_e"
    "Fumarate":   "EX_fum_e",
    "Acetate":    "EX_ac_e",
    "Lactate":    "EX_lac__D_e",
    "Succinate":  "EX_succ_e",
    "Pyruvate":   "EX_pyr_e",
}

def simulate_carbon_source(model, exchange_id, uptake=10.0):
    """Simulate the max growth rate on a single-carbon-source medium.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model (a ``with`` context is used internally,
        so the original model is not modified).
    exchange_id : str
        Exchange reaction id of the carbon source, e.g. "EX_glc__D_e".
    uptake : float
        Max allowed uptake rate (mmol/gDW/h), default 10.

    Returns
    -------
    float
        Max growth rate (1/h) on this carbon source.
    """
    with model:  # with 上下文：退出时自动还原所有改动
        # 关闭所有碳源交换反应的摄取（下界置 0）
        for ex_id in carbon_sources.values():
            model.reactions.get_by_id(ex_id).lower_bound = 0.0
        # 只打开当前碳源，允许最多摄取 uptake
        model.reactions.get_by_id(exchange_id).lower_bound = -uptake
        solution = model.optimize()          # 求解 FBA
        return solution.objective_value      # 返回最大生长速率

print("\n===== Max growth rate on each carbon source "
      "(uptake limit 10 mmol/gDW/h, aerobic) =====")
rows = []
for name, ex_id in carbon_sources.items():
    growth = simulate_carbon_source(model, ex_id)   # 逐个碳源模拟
    rows.append({"carbon_source": name, "exchange_reaction": ex_id,
                 "max_growth_rate_1h": round(growth, 4)})
    print(f"  {name:10s} ({ex_id:15s}) growth rate = {growth:.4f} 1/h")

df_carbon = pd.DataFrame(rows)
carbon_csv = os.path.join(OUTPUT_DIR, "不同碳源生长速率.csv")
df_carbon.to_csv(carbon_csv, index=False, encoding="utf-8-sig")
print("Results saved to:", carbon_csv)

# ----------------------------------------------------------------------
# 2. 葡萄糖摄取速率扫描 —— 简化版表型相平面(PhPP) / 鲁棒性分析
# ----------------------------------------------------------------------
def robustness_scan(model, uptake_values, aerobic=True):
    """Fix a series of glucose uptake rates and compute the max growth rate for each.

    Parameters
    ----------
    model : cobra.Model
        The metabolic model.
    uptake_values : list of float
        Uptake rates to scan (positive values, mmol/gDW/h).
    aerobic : bool
        True keeps the default oxygen supply; False shuts O2 uptake (anaerobic).

    Returns
    -------
    pandas.DataFrame
        Columns: glucose_uptake, max_growth_rate_1h.
    """
    results = []
    with model:                              # 上下文保护，退出自动还原
        if not aerobic:                      # 无氧：关闭 O2 摄取
            model.reactions.EX_o2_e.lower_bound = 0.0
        for uptake in uptake_values:
            # 把下界和上界同时固定为 -uptake，即“强制”按该速率摄取
            model.reactions.EX_glc__D_e.bounds = (-uptake, -uptake)
            sol = model.optimize()
            # 不可行解（如强制摄取过快但无法利用）按 0 处理
            growth = sol.objective_value if sol.status == "optimal" else 0.0
            results.append({"glucose_uptake": round(float(uptake), 2),
                            "max_growth_rate_1h": round(growth, 4)})
    return pd.DataFrame(results)

# 0 到 12 之间均匀取 25 个点进行扫描
uptake_grid = np.linspace(0, 12, 25)

df_aerobic = robustness_scan(model, uptake_grid, aerobic=True)    # 有氧扫描
df_anaerobic = robustness_scan(model, uptake_grid, aerobic=False) # 无氧扫描

print("\n===== Aerobic scan (first 6 rows) =====")
print(df_aerobic.head(6).to_string(index=False))
print("\n===== Anaerobic scan (first 6 rows) =====")
print(df_anaerobic.head(6).to_string(index=False))

# 合并结果并保存 CSV
df_scan = pd.DataFrame({
    "glucose_uptake": uptake_grid.round(2),
    "aerobic_max_growth_1h": df_aerobic["max_growth_rate_1h"],
    "anaerobic_max_growth_1h": df_anaerobic["max_growth_rate_1h"],
})
scan_csv = os.path.join(OUTPUT_DIR, "葡萄糖摄取扫描_表型相平面.csv")
df_scan.to_csv(scan_csv, index=False, encoding="utf-8-sig")
print("\nScan results saved to:", scan_csv)

# ----------------------------------------------------------------------
# 3. 画表型相平面图（生长速率 vs 葡萄糖摄取速率）
# ----------------------------------------------------------------------
fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(df_scan["glucose_uptake"], df_scan["aerobic_max_growth_1h"],
        "o-", label="Aerobic", color="#d62728")
ax.plot(df_scan["glucose_uptake"], df_scan["anaerobic_max_growth_1h"],
        "s--", label="Anaerobic", color="#1f77b4")
ax.set_xlabel("Glucose uptake rate (mmol/gDW/h)")   # x 轴：摄取速率
ax.set_ylabel("Max growth rate (1/h)")              # y 轴：最大生长速率
ax.set_title("Phenotype Phase Plane (robustness analysis)")
ax.legend()
ax.grid(alpha=0.3)
fig_path = os.path.join(OUTPUT_DIR, "表型相平面_葡萄糖摄取扫描.png")
fig.savefig(fig_path, dpi=150, bbox_inches="tight")
print("Figure saved to:", fig_path)

# 解读：有氧曲线斜率更大，说明有氧代谢每摩尔葡萄糖产生更多 ATP、生长更快
print(f"\nAerobic:   growth = {df_scan['aerobic_max_growth_1h'].iloc[-1]} 1/h at uptake 12")
print(f"Anaerobic: growth = {df_scan['anaerobic_max_growth_1h'].iloc[-1]} 1/h at uptake 12")

print("\n=== Script 02_表型预测_生长速率.py finished ===")
