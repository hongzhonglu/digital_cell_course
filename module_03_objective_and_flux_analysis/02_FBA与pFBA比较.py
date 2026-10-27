# -*- coding: utf-8 -*-
"""
模块03 脚本02：FBA 解的简并性 与 FBA/pFBA 比较
================================================
【目的】
1. 理解 FBA 解的"简并性"：达到同一最优目标值的通量分布有无穷多个，
   模型中可能存在"无效循环（futile cycle）"——通量空转但不改变目标值。
2. 用 FVA（通量可变性分析）找出简并反应，构造一个"很差的 FBA 解"。
3. 理解 pFBA（parsimonious FBA，简约通量平衡分析）：
   在保持最优目标值不变的前提下，最小化总通量 sum(|v|)，
   即"用最经济的酶投入实现同样的生长"，剔除无效循环、给出唯一解。

【用法】python 02_FBA与pFBA比较.py
【对应课次】课程5：目标函数与优化原理
【输出】output/02_FBA_vs_pFBA_comparison.png
"""

# ---------------- 导入依赖库 ----------------
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # 无界面环境下保存图片
import matplotlib.pyplot as plt
from cobra.io import load_model, load_json_model, read_sbml_model  # 模型加载
from cobra.flux_analysis import pfba, flux_variability_analysis

# 输出目录（本脚本无绘图，但保留统一约定）
OUTPUT_DIR = os.path.join(os.getcwd(), "module_03_objective_and_flux_analysis/output")
os.makedirs(OUTPUT_DIR, exist_ok=True)
print("Output directory:", OUTPUT_DIR)


# ---------------- 模型加载工具（各模块统一策略，脚本内自带、可独立运行） ----------------
def load_ecoli_core():
    """Load the E. coli core model e_coli_core.

    查找顺序：
      ① 本模块目录下的本地模型文件（e_coli_core.json，BiGG 当前版，
         生物量反应 ID 为 BIOMASS_Ecoli_core_w_GAM）——离线可用，最可靠；
      ② cobra.io.load_model("e_coli_core") 在线下载
         （cobra 0.32 在线下载可能被 BiGG 301 重定向打断）；
      ③ 两者都失败时抛出带明确提示的错误。
    """
    __file__ ='module_03_objective_and_flux_analysis'
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


def find_biomass_rxn(model):
    """Auto-detect the biomass reaction, compatible with both namings:
    'Biomass_Ecoli_core' (cobra textbook version) and
    'BIOMASS_Ecoli_core_w_GAM' (current BiGG version)."""
    for r in model.reactions:
        if "biomass" in r.id.lower():      # id 小写后含 biomass 即为生物量反应
            return r
    # 兜底：取当前目标函数中出现的第一个反应
    for r in model.reactions:
        if r.objective_coefficient != 0:
            return r
    raise ValueError("No biomass reaction found in the model; check the model file.")



# ---------------- 加载模型 ----------------
model = load_ecoli_core()
BIOMASS = find_biomass_rxn(model).id  # 自动探测生物量反应 ID（两套命名通用）

# ============================================================
# 1. 基础 FBA：最优生长率
# ============================================================
fba_sol = model.optimize()
print("=" * 60)
print(f"[1] Wild-type FBA optimal growth rate = {fba_sol.objective_value:.4f} 1/h")
print(f"    Total flux sum(|v|) of this solution = {fba_sol.fluxes.abs().sum():.2f}")


# ============================================================
# 2. 构造一个"很差的 FBA 解"：强制循环高通量，生长率不变
# ============================================================
print("=" * 60)
print("[3] Force SUCDi flux = 400 (spinning cycle) and re-run FBA:")
with model:  # with 语句退出时自动还原模型，适合课堂演示
    model.reactions.SUCDi.lower_bound = 400  # 强制无效循环携带 400 的通量
    bad_fba = model.optimize()
    print(f"    Growth rate is still {bad_fba.objective_value:.4f} 1/h (objective unchanged!)")
    print(f"    But total flux sum(|v|) inflates to {bad_fba.fluxes.abs().sum():.2f}")
    print("    Biological meaning: the cell spends lots of 'enzyme' spinning a cycle")
    print("    with zero growth benefit.")

# ============================================================
# 3. pFBA：同一最优目标下的最经济解
# ============================================================
print("=" * 60)
print("[4] pFBA on the 'forced-cycle' model:")
with model:
    model.reactions.SUCDi.lower_bound = 0  # 只要求生长最优，不强制循环
    pfba_sol = pfba(model)
    growth_pfba = pfba_sol.fluxes[BIOMASS]
    print(f"    pFBA growth rate = {growth_pfba:.4f} 1/h (same as FBA)")
    print(f"    pFBA total flux sum(|v|) = {pfba_sol.fluxes.abs().sum():.2f}")
    print(f"    SUCDi flux in pFBA = {pfba_sol.fluxes['SUCDi']:.4f}"
          " (back to the physiological minimum)")

# ============================================================
# 4. 比较三条解：普通FBA / 强制循环FBA / pFBA
# ============================================================
labels = ["FBA (plain)", "FBA (forced cycle)", "pFBA"]
growth = [fba_sol.objective_value, bad_fba.objective_value, growth_pfba]
totals = [fba_sol.fluxes.abs().sum(), bad_fba.fluxes.abs().sum(),
          pfba_sol.fluxes.abs().sum()]

print("=" * 60)
print("[5] Summary of the three solutions:")
print(pd.DataFrame({"growth rate": [round(g, 4) for g in growth],
                    "sum(|v|)": [round(t, 2) for t in totals]}, index=labels))

# 挑几条代表性反应比较通量差异（含循环反应 SUCDi/FRD7 与核心反应）
show_ids = ["SUCDi", "FRD7", "MDH", "ATPS4r", "EX_glc__D_e", "EX_o2_e"]
cmp_df = pd.DataFrame({
    "FBA (forced cycle)": bad_fba.fluxes[show_ids].round(3),
    "pFBA": pfba_sol.fluxes[show_ids].round(3),
})
print("\nFlux comparison for representative reactions:")
print(cmp_df)

# ============================================================
# 5. 画对比条形图（英文标注）
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# 左图：总通量对比（生长率都相同）
colors = ["#7f7f7f", "#d62728", "#1f77b4"]
axes[0].bar(labels, totals, color=colors)
axes[0].set_ylabel("sum(|v|)  (mmol/gDW/h)")
axes[0].set_title("Total flux (all three have the same growth rate = 0.874 1/h)")
for i, v in enumerate(totals):
    axes[0].text(i, v, f"{v:.0f}", ha="center", va="bottom")

# 右图：代表反应通量对比（强制循环FBA vs pFBA）
x = range(len(show_ids))
width = 0.38
axes[1].bar([i - width / 2 for i in x], bad_fba.fluxes[show_ids],
            width=width, label="FBA (forced cycle)", color="#d62728")
axes[1].bar([i + width / 2 for i in x], pfba_sol.fluxes[show_ids],
            width=width, label="pFBA", color="#1f77b4")
axes[1].set_xticks(list(x))
axes[1].set_xticklabels(show_ids, rotation=30, ha="right")
axes[1].set_ylabel("Flux (mmol/gDW/h)")
axes[1].set_title("Reaction-by-reaction: pFBA removes the futile cycle")
axes[1].legend()

plt.tight_layout()
png_path = os.path.join(OUTPUT_DIR, "02_FBA_vs_pFBA_comparison.png")
plt.savefig(png_path, dpi=150)
print("=" * 60)
print("Figure saved:", png_path)
