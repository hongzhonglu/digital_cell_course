# -*- coding: utf-8 -*-
"""
模块04 脚本02：通量上下界与方向约束
================================================
【目的】
1. 理解反应边界 (lower_bound, upper_bound) 是 GEM 最直接的约束语言：
   - lb < 0, ub > 0：可逆反应；
   - lb = 0：不可逆反应（方向约束来自热力学/反应机制）；
   - lb = ub = 0：反应关闭（等价于敲除）。
2. 演示收紧/放开关键反应边界对最优解的影响（可行域↔目标值）。
3. 重点演示 ATPM（ATP 维持消耗）下界的意义：
   它强制细胞分出部分资源用于"非生长维持"，下界越大，可生长越少。
4. 演示"收紧不可逆性约束"（如强制某可逆反应单向）会如何改变解。

【用法】python 02_通量上下界与方向约束.py
【对应课次】课程6：GEMs运行所需的约束
【输出】output/02_bound_constraint_effects.png
"""

# ---------------- 导入依赖库 ----------------
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cobra.io import load_model, load_json_model, read_sbml_model  # 模型加载

# 图输出目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)


# ---------------- 模型加载工具（各模块统一策略，脚本内自带、可独立运行） ----------------
def load_ecoli_core():
    """Load e_coli_core: (1) local file in this module directory
    (e_coli_core.json/xml, current BiGG version whose biomass reaction is
    BIOMASS_Ecoli_core_w_GAM) -> (2) cobra online download (may be broken by
    a BiGG 301 redirect) -> (3) raise a clear error."""
    here = os.path.dirname(os.path.abspath(__file__))    # 本脚本所在目录
    local_json = os.path.join(here, "e_coli_core.json")  # 本地 JSON 模型
    local_xml = os.path.join(here, "e_coli_core.xml")    # 本地 SBML 模型（备选）
    if os.path.exists(local_json):
        return load_json_model(local_json)               # ① 优先本地文件
    if os.path.exists(local_xml):
        return read_sbml_model(local_xml)                # ① 优先本地文件
    try:
        return load_model("e_coli_core")                 # ② 在线下载兜底
    except Exception as exc:
        raise RuntimeError(                              # ③ 明确报错
            "Failed to load e_coli_core: local e_coli_core.json is missing "
            "and online download failed. Copy e_coli_core.json into this "
            "module directory and retry."
        ) from exc


# ---------------- 加载模型 ----------------
model = load_ecoli_core()
wt_growth = model.optimize().objective_value
print("=" * 60)
print(f"Wild-type (default bounds) growth rate = {wt_growth:.4f} 1/h")

# ============================================================
# 1. 查看典型反应的边界：可逆 / 不可逆 / 固定
# ============================================================
print("=" * 60)
print("[1] Bounds of typical reactions:")
examples = ["PGI", "FBA", "ATPM", "EX_glc__D_e", "MDH", "SUCDi"]
for rid in examples:
    r = model.reactions.get_by_id(rid)
    direction = "reversible" if r.lower_bound < 0 < r.upper_bound else \
                ("irreversible (forward)" if r.lower_bound >= 0
                 else "irreversible (reverse)")
    print(f"    {rid:<12s} bounds = {str(r.bounds):<18s} {direction}")
print("    Reading: direction constraints (irreversibility) come from reaction")
print("    thermodynamics - an important prior constraint of a GEM.")

# ============================================================
# 2. ATPM 下界：非生长维持能耗约束
# ============================================================
print("=" * 60)
print("[2] Scanning the ATPM (ATP maintenance) lower bound:")
print(f"    ATPM default bounds = {model.reactions.ATPM.bounds} (lower bound 8.39 =")
print("    measured non-growth-associated maintenance cost of E. coli)")
atpm_lbs = [0, 8.39, 15, 25, 35, 45, 55]
atpm_growth = []
for lb in atpm_lbs:
    with model:
        model.reactions.ATPM.lower_bound = lb
        g = model.optimize().objective_value
    atpm_growth.append(g)
    print(f"    ATPM lb = {lb:>6.2f}  ->  growth = {g:.4f} 1/h")
print("    -> The higher the maintenance cost, the less ATP is left for")
print("       biomass synthesis, so growth falls monotonically;")
print("       once the bound makes the LP infeasible, the model predicts death.")

# 找 ATPM 致死临界值（LP 不可行：所有 ATP 都用于维持也供不起）
print("\n    Probing the lethal threshold (theoretical max ~175 in this model):")
for lb in [100, 150, 170, 180]:
    with model:
        model.reactions.ATPM.lower_bound = lb
        s = model.optimize()
    print(f"    ATPM lb = {lb:>3}  ->  status = {s.status:>9s}, growth = {s.objective_value}")

# ============================================================
# 3. 收紧关键反应上界：以呼吸链 CYTBD 为例
# ============================================================
print("=" * 60)
print("[3] Tightening the CYTBD upper bound (limiting respiratory capacity):")
cytbd_ubs = [1000, 40, 30, 20, 10, 0]
cytbd_growth, cytbd_ac = [], []
for ub in cytbd_ubs:
    with model:
        model.reactions.CYTBD.upper_bound = ub
        s = model.optimize()
    cytbd_growth.append(s.objective_value)
    cytbd_ac.append(s.fluxes["EX_ac_e"])
    print(f"    CYTBD ub = {ub:>6.1f}  ->  growth = {s.objective_value:.4f},"
          f" acetate secretion = {s.fluxes['EX_ac_e']:.3f}")
print("    -> With respiration limited the cell shifts to fermentation")
print("       (acetate overflow), similar to overflow metabolism.")

# ============================================================
# 4. 收紧方向约束：把可逆反应 PGI 强制为单向
# ============================================================
print("=" * 60)
print("[4] Direction-constraint demo:")
print("    PGI (phosphoglucose isomerase) is reversible by default.")
print("    Forcing it forward-only (lb=0) mimics a tighter thermodynamic")
print("    constraint - the feasible region shrinks, so the optimum can only")
print("    stay the same or get worse, never better:")
with model:
    model.reactions.PGI.lower_bound = 0  # 强制 G6P→F6P 单向
    s = model.optimize()
print(f"    Growth with PGI forced forward-only = {s.objective_value:.4f} 1/h"
      f" (default {wt_growth:.4f})")
print("    -> Tightening a direction constraint shrinks the feasible region;")
print("       the optimal objective stays unchanged or worsens.")

# ============================================================
# 5. 汇总绘图（英文标注）
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

axes[0].plot(atpm_lbs, atpm_growth, "o-", color="#d62728")
axes[0].axvline(8.39, ls="--", color="gray", label="default 8.39")
axes[0].set_xlabel("ATPM lower bound (mmol/gDW/h)")
axes[0].set_ylabel("Max growth rate (1/h)")
axes[0].set_title("ATP maintenance constraint vs growth rate")
axes[0].legend()

axes[1].plot(cytbd_ubs, cytbd_growth, "s-", color="#1f77b4", label="growth rate")
ax2 = axes[1].twinx()
ax2.plot(cytbd_ubs, cytbd_ac, "^--", color="#ff7f0e", label="acetate secretion")
axes[1].set_xlabel("CYTBD upper bound (mmol/gDW/h)")
axes[1].set_ylabel("Max growth rate (1/h)", color="#1f77b4")
ax2.set_ylabel("Acetate secretion (mmol/gDW/h)", color="#ff7f0e")
axes[1].set_title("Respiratory capacity constraint -> growth and overflow")
lines = axes[1].get_lines() + ax2.get_lines()
axes[1].legend(lines, [l.get_label() for l in lines], loc="center right")

plt.tight_layout()
png_path = os.path.join(OUTPUT_DIR, "02_bound_constraint_effects.png")
plt.savefig(png_path, dpi=150)
print("=" * 60)
print("Figure saved:", png_path)
print("[Summary] Bounds = constraints: tightening an upper/lower bound only")
print("  shrinks the feasible region and keeps or lowers the optimum;")
print("  the ATPM lower bound encodes 'maintenance energy'; direction")
print("  constraints encode 'thermodynamic irreversibility'.")
print("Script 02 finished.")
