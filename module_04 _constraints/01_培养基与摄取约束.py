# -*- coding: utf-8 -*-
"""
模块04 脚本01：培养基与摄取约束
================================================
【目的】
1. 理解 GEM 预测的前提是"环境约束"：培养基由交换反应（EX_*）的下界定义，
   model.medium 给出当前各营养物的最大摄取速率。
2. 学会逐个改变葡萄糖 / 氧气 / 氨的摄取下界，观察最大生长率的变化，
   体会"约束收紧 → 可行域缩小 → 最优目标不增"的基本规律。
3. 模拟厌氧条件（EX_o2_e 下界=0），比较有氧/厌氧下的产物分泌谱
   （乙醇、乙酸、甲酸、琥珀酸等），理解氧化还原约束对发酵产物分布的决定作用。

【用法】python 01_培养基与摄取约束.py
【对应课次】课程6：GEMs运行所需的约束
【输出】output/01_medium_and_secretion_comparison.png
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


# ---------------- 加载模型并查看默认培养基 ----------------
model = load_ecoli_core()
print("=" * 60)
print("[1] Default medium, model.medium (exchange ID: max uptake rate):")
print(model.medium)
print("Reading: EX_glc__D_e = 10 means glucose uptake is capped at 10 mmol/gDW/h;")
print("         EX_o2_e = 1000 means O2 is practically unlimited (internal")
print("         constraints do the real limiting).")
# medium 与反应下界的关系：medium 中的值会写入对应交换反应的 lower_bound
print("\nEX_glc__D_e bounds:", model.reactions.EX_glc__D_e.bounds)
print("EX_o2_e bounds:", model.reactions.EX_o2_e.bounds)
print("EX_nh4_e bounds:", model.reactions.EX_nh4_e.bounds)

wt_growth = model.optimize().objective_value
print(f"\nMax growth rate on the default medium = {wt_growth:.4f} 1/h")

# ============================================================
# 2. 逐个改变葡萄糖摄取下界（碳源限制）
# ============================================================
print("=" * 60)
print("[2] Glucose uptake scan (EX_glc__D_e lower bound = -10 -> 0):")
glc_uptakes = [10, 8, 6, 4, 2, 1, 0.5, 0]
glc_growth = []
for u in glc_uptakes:
    with model:  # with 块内修改边界，退出自动还原
        model.reactions.EX_glc__D_e.lower_bound = -u
        g = model.optimize().objective_value
        # 无碳源时 ATPM 强制维持消耗会使"生长"出现微小负值，
        # 生物学上解读为 0（细胞不生长），这里截断显示
        g = max(g, 0.0)
    glc_growth.append(g)
    print(f"    Glucose uptake <= {u:>5.1f}  ->  growth = {g:.4f} 1/h")
print("    -> Growth drops linearly to 0 with carbon uptake: carbon is the")
print("       rate-limiting substrate for growth.")

# ============================================================
# 3. 改变氧气摄取下界（从有氧到厌氧）
# ============================================================
print("=" * 60)
print("[3] O2 uptake scan (EX_o2_e lower bound = -1000 -> 0):")
o2_uptakes = [1000, 20, 15, 10, 5, 2, 0]
o2_growth = []
for u in o2_uptakes:
    with model:
        model.reactions.EX_o2_e.lower_bound = -u
        g = model.optimize().objective_value
    o2_growth.append(g)
    print(f"    O2 uptake <= {u:>6.1f}  ->  growth = {g:.4f} 1/h")
print("    -> Less O2 means less respiratory ATP, so growth drops;")
print("       anaerobic growth is still possible (fermentation).")

# ============================================================
# 4. 改变氨摄取下界（氮源限制）
# ============================================================
print("=" * 60)
print("[4] NH4 uptake scan (EX_nh4_e lower bound = -1000 -> 0):")
nh4_uptakes = [1000, 10, 5, 4, 3, 2, 1, 0]
nh4_growth = []
for u in nh4_uptakes:
    with model:
        model.reactions.EX_nh4_e.lower_bound = -u
        g = model.optimize().objective_value
    nh4_growth.append(g)
    print(f"    NH4 uptake <= {u:>6.1f}  ->  growth = {g:.4f} 1/h")
print("    -> Nitrogen limitation restricts growth too (protein synthesis")
print("       needs nitrogen).")

# ============================================================
# 5. 厌氧条件模拟：比较有氧 vs 厌氧的产物分泌谱
# ============================================================
print("=" * 60)
print("[5] Aerobic vs anaerobic: secretion fluxes of major fermentation products:")
secretion_ids = ["EX_ac_e", "EX_etoh_e", "EX_for_e", "EX_succ_e",
                 "EX_lac__D_e", "EX_co2_e"]
# 交换反应 ID 对应的产物英文名（用于表格与绘图标注）
product_name = {"EX_ac_e": "acetate", "EX_etoh_e": "ethanol",
                "EX_for_e": "formate", "EX_succ_e": "succinate",
                "EX_lac__D_e": "lactate", "EX_co2_e": "CO2"}

# --- 有氧解 ---
sol_aerobic = model.optimize()
# --- 厌氧解：氧气摄取下界置 0 ---
with model:
    model.reactions.EX_o2_e.lower_bound = 0
    sol_anaerobic = model.optimize()

cmp = pd.DataFrame({
    "product": [product_name[i] for i in secretion_ids],
    "aerobic": sol_aerobic.fluxes[secretion_ids].round(3).values,
    "anaerobic": sol_anaerobic.fluxes[secretion_ids].round(3).values,
}, index=secretion_ids)
print(cmp)
print(f"\nAerobic growth = {sol_aerobic.objective_value:.4f} 1/h, "
      f"anaerobic growth = {sol_anaerobic.objective_value:.4f} 1/h")
print("-> Without O2 the respiratory chain cannot drain electrons; the cell")
print("   regenerates NAD+ by secreting reduced products (ethanol/acetate),")
print("   and CO2 secretion changes markedly - the redox constraint at work.")

# ============================================================
# 6. 可视化（英文标注）
# ============================================================
fig, axes = plt.subplots(1, 3, figsize=(16, 5))

# 左图：三种底物摄取-生长曲线（统一x为摄取上限）
axes[0].plot(glc_uptakes, glc_growth, "o-", label="glucose", color="#d62728")
axes[0].plot(o2_uptakes, o2_growth, "s-", label="O2", color="#1f77b4")
axes[0].plot(nh4_uptakes, nh4_growth, "^-", label="NH4", color="#2ca02c")
axes[0].set_xlabel("Uptake upper bound (mmol/gDW/h)")
axes[0].set_ylabel("Max growth rate (1/h)")
axes[0].set_title("Substrate uptake constraint -> growth rate")
axes[0].legend()
axes[0].set_xscale("symlog", linthresh=1)

# 中图：葡萄糖限制下生长率（线性关系特写）
axes[1].plot(glc_uptakes, glc_growth, "o-", color="#d62728")
axes[1].set_xlabel("Glucose uptake upper bound (mmol/gDW/h)")
axes[1].set_ylabel("Max growth rate (1/h)")
axes[1].set_title("Carbon limitation: ~linear response")

# 右图：有氧/厌氧分泌对比
x = range(len(secretion_ids))
width = 0.38
axes[2].bar([i - width / 2 for i in x], sol_aerobic.fluxes[secretion_ids],
            width=width, label=f"aerobic (mu={sol_aerobic.objective_value:.2f})",
            color="#1f77b4")
axes[2].bar([i + width / 2 for i in x], sol_anaerobic.fluxes[secretion_ids],
            width=width, label=f"anaerobic (mu={sol_anaerobic.objective_value:.2f})",
            color="#ff7f0e")
axes[2].set_xticks(list(x))
axes[2].set_xticklabels([product_name[i] for i in secretion_ids],
                        rotation=30, ha="right")
axes[2].set_ylabel("Secretion flux (mmol/gDW/h)")
axes[2].set_title("Aerobic vs anaerobic: fermentation product spectrum")
axes[2].legend()

plt.tight_layout()
png_path = os.path.join(OUTPUT_DIR, "01_medium_and_secretion_comparison.png")
plt.savefig(png_path, dpi=150)
print("=" * 60)
print("Figure saved:", png_path)
print("[Summary] Medium = lower bounds of exchange reactions; tightening")
print("  substrate constraints only shrinks the feasible region and lowers")
print("  growth; the O2 constraint reroutes electrons and shapes the")
print("  fermentation product spectrum.")
print("Script 01 finished.")
