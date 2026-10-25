# -*- coding: utf-8 -*-
"""
模块03 脚本04：通量可变性分析（FVA）——野生型 vs pgi 敲除株
================================================
【目的】
1. 掌握 cobra.flux_analysis.flux_variability_analysis（FVA）的用法：
   对每个反应分别求"在目标函数保持一定比例最优值的前提下"的最小/最大通量，
   得到每条反应的通量可行区间 [min, max]。
2. 理解 FVA 与 FBA/pFBA 回答的是"不同的问题"：
   - FBA  ：给出"一个"最优目标值和"一个"最优通量分布（解通常不唯一）；
   - pFBA ：在所有最优解中挑出总通量最小的"一个"代表解（仍是一个点解）；
   - FVA  ：不挑代表解，而是刻画"整个最优解空间"——回答
            "某反应的通量在不牺牲目标值的情况下能在多大范围内变化？"
            即解的"唯一性/柔性（flexibility）"问题。
   因此 FVA 是衡量 FBA 预测置信度的工具：区间越窄，预测越确定。
3. 实操：野生型 vs pgi(b4025) 敲除株（fraction_of_optimum=0.9，即允许
   牺牲至多 10% 生长率换取通量柔性），对比关键反应的通量区间；
   统计"关闭（区间≈0）/ 固定（区间宽度≈0）/ 可变"反应数量；
   绘制通量区间对比图（英文标注）保存为 png。

【用法】python 04_flux_variability_analysis.py
【对应课次】课程5：目标函数与优化原理
【输出】
  output/04_FVA_wildtype_vs_dpgi.png   关键反应通量区间对比图
  output/04_fva_wildtype.csv           野生型 FVA 结果表
  output/04_fva_dpgi.csv               Δpgi 敲除株 FVA 结果表
"""

# ---------------- 导入依赖库 ----------------
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # 无界面环境下保存图片
import matplotlib.pyplot as plt
from cobra.io import load_model, load_json_model, read_sbml_model  # 模型加载
from cobra.flux_analysis import flux_variability_analysis

# 图输出目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# FVA 中允许牺牲的目标比例：0.9 = 保留至少 90% 最优生长率
# 设为 0.9 而非 1.0 的原因：允许少量次优，能看到敲除株"柔性"如何变化，
# 教学上更能体现"解空间"的概念
FRACTION = 0.9
# 判定"区间宽度≈0"的数值容差（LP 求解存在数值误差，不能直接用 == 0）
ZERO_TOL = 1e-6


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


def find_biomass_rxn(model):
    """Auto-detect the biomass reaction, compatible with 'Biomass_Ecoli_core'
    (textbook version) and 'BIOMASS_Ecoli_core_w_GAM' (current BiGG version)."""
    for r in model.reactions:
        if "biomass" in r.id.lower():      # id 小写后含 biomass 即为生物量反应
            return r
    # 兜底：取当前目标函数中出现的第一个反应
    for r in model.reactions:
        if r.objective_coefficient != 0:
            return r
    raise ValueError("No biomass reaction found in the model; check the model file.")


def classify_reactions(fva_df, tol=ZERO_TOL):
    """Classify reactions from an FVA result table into three groups.

    分类规则：
      - blocked（关闭）：min≈max≈0，该反应在可行解空间中根本无法携带通量；
      - fixed（固定）：min≈max≠0，通量被约束/目标唯一确定（预测最可信）；
      - variable（可变）：min 与 max 有实质差距（解不唯一，存在柔性）。

    Parameters
    ----------
    fva_df : pandas.DataFrame
        FVA 结果，含 'minimum' 与 'maximum' 两列，index 为反应 ID。
    tol : float
        判定区间宽度/通量值为 0 的数值容差。

    Returns
    -------
    pandas.Series
        每个反应的分类标签：'blocked' / 'fixed' / 'variable'。
    """
    width = fva_df["maximum"] - fva_df["minimum"]        # 区间宽度
    magnitude = fva_df[["minimum", "maximum"]].abs().max(axis=1)  # 区间最大绝对值
    labels = pd.Series("variable", index=fva_df.index)   # 默认可变
    labels[magnitude <= tol] = "blocked"                 # 整个区间都在 0 附近
    labels[(magnitude > tol) & (width <= tol)] = "fixed"  # 非零但宽度为 0
    return labels


# ---------------- 加载模型 ----------------
model = load_ecoli_core()
BIOMASS = find_biomass_rxn(model).id  # 自动探测生物量反应 ID（两套命名通用）
wt_growth = model.optimize().objective_value
print("=" * 60)
print(f"Wild-type optimal growth rate = {wt_growth:.4f} 1/h")
print(f"FVA setting: fraction_of_optimum = {FRACTION}"
      f" (growth is allowed to drop to {FRACTION*100:.0f}% of optimum)")

# ============================================================
# 1. 野生型 FVA
# ============================================================
# FVA 的本质：对每个反应 v_i 依次解两个 LP：
#   min v_i 和 max v_i，约束为 S·v=0、lb<=v<=ub、目标值 >= fraction*最优值
# 95 个反应 × 2 个 LP = 190 次求解（核心模型秒级完成）
print("=" * 60)
print("[1] Running FVA on the wild-type model (190 LPs) ...")
fva_wt = flux_variability_analysis(model, fraction_of_optimum=FRACTION)
print("    Done. Preview (first 5 rows):")
print(fva_wt.head().round(4))

# 统计三类反应数量：blocked / fixed / variable
cls_wt = classify_reactions(fva_wt)
print("\n    Reaction classes in wild type:")
for cls in ["blocked", "fixed", "variable"]:
    print(f"      {cls:<9s}: {(cls_wt == cls).sum():>3d} reactions")

# ============================================================
# 2. pgi 敲除株 FVA
# ============================================================
# 用 with model: 上下文做敲除，退出后自动还原模型，不污染野生型
print("=" * 60)
print("[2] Knocking out pgi (b4025) and running FVA on the mutant ...")
with model:
    model.genes.b4025.knock_out()          # 按 GPR 规则关闭 PGI 反应
    ko_growth = model.optimize().objective_value
    print(f"    dpgi optimal growth rate = {ko_growth:.4f} 1/h"
          f" ({ko_growth/wt_growth*100:.1f}% of wild type)")
    fva_ko = flux_variability_analysis(model, fraction_of_optimum=FRACTION)
    cls_ko = classify_reactions(fva_ko)
    print("\n    Reaction classes in dpgi:")
    for cls in ["blocked", "fixed", "variable"]:
        print(f"      {cls:<9s}: {(cls_ko == cls).sum():>3d} reactions")

# 解读：敲除 pgi 后 G6P 只能走戊糖磷酸途径，
# 原本"可变"的一些反应被新约束"钉死"或彻底失去通量（→ blocked），
# 解空间整体变窄——这就是"柔性下降"
n_to_fixed = ((cls_wt == "variable") & (cls_ko == "fixed")).sum()
n_to_blocked = ((cls_wt != "blocked") & (cls_ko == "blocked")).sum()
print(f"\n    dpgi effect: {n_to_fixed} reactions changed variable -> fixed,"
      f" {n_to_blocked} reactions became blocked.")

# 柔性变化的量化：比较两株中每个反应的区间宽度之和与"显著收窄"的反应数
# 区间宽度 width = max - min；宽度收窄即该反应的通量"不确定性"下降
width_wt = fva_wt["maximum"] - fva_wt["minimum"]
width_ko = fva_ko["maximum"] - fva_ko["minimum"]
n_shrunk = ((width_wt - width_ko) > 1).sum()   # 宽度收窄超过 1 的反应数
n_grown = ((width_ko - width_wt) > 1).sum()    # 宽度变宽超过 1 的反应数
print(f"    Flux-range widths: total {width_wt.sum():.1f} (WT) ->"
      f" {width_ko.sum():.1f} (dpgi); {n_shrunk} reactions narrowed by >1,"
      f" {n_grown} widened.")
print("    Interpretation: losing PGI forces carbon through the pentose")
print("    phosphate pathway; PPP/glycolysis ranges (PGI, G6PDH2r, GND, RPE...)")
print("    narrow the most - the solution space shrinks (less flexibility).")

# 保存两个 FVA 结果表，供课后复查
csv_wt = os.path.join(OUTPUT_DIR, "04_fva_wildtype.csv")
csv_ko = os.path.join(OUTPUT_DIR, "04_fva_dpgi.csv")
fva_wt.assign(category=cls_wt).round(6).to_csv(csv_wt, index_label="reaction")
fva_ko.assign(category=cls_ko).round(6).to_csv(csv_ko, index_label="reaction")
print("    FVA tables saved:", csv_wt, "and", csv_ko)

# ============================================================
# 3. 关键反应的通量区间对比
# ============================================================
# 挑一组能体现"糖酵解 vs 戊糖磷酸途径改道"的代表反应：
#   PGI（被敲除，应 blocked）、G6PDH2r（PPP 入口，敲除后通量被迫增大）、
#   PFK（糖酵解下游）、Biomass（目标）、EX_ac_e / EX_co2_e（分泌）
key_ids = [BIOMASS, "PGI", "G6PDH2r", "PFK", "GAPD", "TKT1",
           "EX_glc__D_e", "EX_ac_e", "EX_co2_e", "EX_o2_e"]
key_ids = [r for r in key_ids if r in fva_wt.index]  # 防御：只保留存在的反应

cmp = pd.DataFrame({
    "WT min": fva_wt.loc[key_ids, "minimum"],
    "WT max": fva_wt.loc[key_ids, "maximum"],
    "dpgi min": fva_ko.loc[key_ids, "minimum"],
    "dpgi max": fva_ko.loc[key_ids, "maximum"],
}).round(3)
print("=" * 60)
print("[3] Flux ranges of key reactions (wild type vs dpgi):")
print(cmp)

# ============================================================
# 4. 画通量区间对比图（英文标注）
# ============================================================
# 每个反应画两条水平"区间棒"：上=野生型（蓝），下=Δpgi（橙），
# 棒的长度即 FVA 区间 [min, max]，直观对比两株的柔性差异
fig, ax = plt.subplots(figsize=(11, 6))
y_pos = list(range(len(key_ids)))[::-1]   # 反应从上到下排列
bar_h = 0.32                              # 区间棒的半高偏移

for y, rid in zip(y_pos, key_ids):
    # 野生型区间棒
    ax.plot([fva_wt.loc[rid, "minimum"], fva_wt.loc[rid, "maximum"]],
            [y + bar_h / 2, y + bar_h / 2],
            color="#1f77b4", lw=8, solid_capstyle="butt",
            label="Wild type" if rid == key_ids[0] else None)
    # Δpgi 区间棒
    ax.plot([fva_ko.loc[rid, "minimum"], fva_ko.loc[rid, "maximum"]],
            [y - bar_h / 2, y - bar_h / 2],
            color="#ff7f0e", lw=8, solid_capstyle="butt",
            label="dpgi" if rid == key_ids[0] else None)

ax.axvline(0, color="gray", ls="--", lw=1)  # 零通量参考线
ax.set_yticks(y_pos)
ax.set_yticklabels(key_ids)
ax.set_xlabel("Flux range (mmol/gDW/h)")
ax.set_title(f"FVA flux ranges: wild type vs dpgi "
             f"(fraction_of_optimum = {FRACTION})")
ax.legend(loc="lower right")
plt.tight_layout()
png_path = os.path.join(OUTPUT_DIR, "04_FVA_wildtype_vs_dpgi.png")
plt.savefig(png_path, dpi=150)
print("Figure saved:", png_path)

# ============================================================
# 5. 小结
# ============================================================
print("=" * 60)
print("[Summary] FVA vs FBA/pFBA: different questions")
print("  FBA  -> one optimal value + one (non-unique) flux distribution;")
print("  pFBA -> one unique, most 'economical' optimal distribution;")
print("  FVA  -> the FULL range of each reaction across (near-)optimal")
print("          solutions, i.e. solution uniqueness / flexibility.")
print("  Wide range = the model cannot pin down that flux (low confidence);")
print("  zero-width range = the prediction is uniquely determined.")
print("  After the pgi knockout the solution space shrinks: fewer variable")
print("  reactions, more fixed/blocked ones (less metabolic flexibility).")
print("Script 04 finished.")
