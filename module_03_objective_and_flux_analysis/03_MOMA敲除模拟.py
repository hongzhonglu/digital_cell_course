# -*- coding: utf-8 -*-
"""
模块03 脚本03：MOMA 基因敲除模拟
================================================
【目的】
1. 理解三种预测敲除株表型的假设差异：
   - FBA ：假设敲除株"瞬间进化"到新的生长最优（最大化生物量）。
           往往高估敲除株生长率。
   - pFBA：在 FBA 最优基础上最小化总通量，得到经济的野生型参考解。
   - MOMA（Minimization Of Metabolic Adjustment，代谢调整最小化）：
           假设敲除株尚未进化，其通量分布尽可能"接近野生型"，
           即最小化 ||v_mutant - v_wildtype||，常用 L1（线性）或 L2（二次）范数。
           更贴近实验上"刚敲除、未驯化"的菌株。
2. 实操：以野生型 pFBA 解为参考，对 gltA(b0720) 与 pgi(b4025) 做敲除，
   比较 FBA 敲除预测 vs MOMA 敲除预测的生长率。
3. 说明 quadratic MOMA 需要 QP 求解器（如 Gurobi/CPLEX），
   本环境无 QP 求解器，故采用 linear MOMA（结论方向一致）。

【用法】python 03_MOMA敲除模拟.py
【对应课次】课程5：目标函数与优化原理
【输出】output/03_MOMA_knockout_comparison.png
"""

# ---------------- 导入依赖库 ----------------
import os
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from cobra.io import load_model, load_json_model, read_sbml_model  # 模型加载
from cobra.flux_analysis import pfba, moma

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



# ---------------- 加载模型、求野生型参考解 ----------------
model = load_ecoli_core()
BIOMASS = find_biomass_rxn(model).id  # 自动探测生物量反应 ID（两套命名通用）
wt_fba = model.optimize()
wt_pfba = pfba(model)  # MOMA 的标准做法：以野生型 pFBA 解作为参考通量分布
wt_growth = wt_fba.objective_value
print("=" * 60)
print(f"Wild-type FBA growth rate = {wt_growth:.4f} 1/h")
print(f"Wild-type pFBA growth rate = {wt_pfba.fluxes[BIOMASS]:.4f} 1/h"
      " (used as the MOMA reference solution)")


# ============================================================
# 定义敲除模拟函数
# ============================================================
def simulate_knockout(model, gene_id, reference):
    """Knock out a gene and predict mutant growth with FBA and linear MOMA.

    Parameters
    ----------
    model : cobra.Model
        Wild-type model (copied internally, so the original is untouched).
    gene_id : str
        ID of the gene to knock out (e.g. 'b4025').
    reference : cobra.Solution
        Wild-type reference flux distribution (usually the pFBA solution).

    Returns
    -------
    dict with gene name, FBA-predicted growth, MOMA-predicted growth,
    and the two solution objects.
    """
    ko = model.copy()                        # 复制模型，避免修改野生型
    gene = ko.genes.get_by_id(gene_id)
    reactions = [r.id for r in gene.reactions]
    gene.knock_out()                         # 依据 GPR 规则关闭相关反应
    # --- 方法1：FBA 假设（敲除株重新进化到生长最优）---
    fba_ko = ko.optimize()
    # --- 方法2：linear MOMA 假设（通量尽量接近野生型参考解）---
    moma_ko = moma(ko, solution=reference, linear=True)
    return {
        "gene_id": gene_id,
        "gene_name": gene.name,
        "affected_reactions": ",".join(reactions),
        "fba_growth": fba_ko.objective_value,
        "moma_growth": moma_ko.fluxes[BIOMASS],
        "_fba_sol": fba_ko, "_moma_sol": moma_ko,
    }


# ============================================================
# 1. 敲除 gltA (b0720)：TCA 循环柠檬酸合酶
# ============================================================
print("=" * 60)
print("[1] Knocking out gltA (b0720), which catalyzes CS (citrate synthase):")
res1 = simulate_knockout(model, "b0720", wt_pfba)
print(f"    FBA  predicted growth = {res1['fba_growth']:.4f} 1/h")
print(f"    MOMA predicted growth = {res1['moma_growth']:.4f} 1/h")
print("    -> Both are ~0: gltA is essential under this condition;")
print("       no assumption can rescue it (no bypass for alpha-ketoglutarate).")

# ============================================================
# 2. 敲除 pgi (b4025)：糖酵解磷酸葡萄糖异构酶——经典 MOMA 案例
# ============================================================
print("=" * 60)
print("[2] Knocking out pgi (b4025), which catalyzes PGI (phosphoglucose isomerase):")
res2 = simulate_knockout(model, "b4025", wt_pfba)
print(f"    Wild type             = {wt_growth:.4f} 1/h")
print(f"    FBA  predicted growth = {res2['fba_growth']:.4f} 1/h"
      f" ({res2['fba_growth']/wt_growth*100:.1f}% retained)")
print(f"    MOMA predicted growth = {res2['moma_growth']:.4f} 1/h"
      f" ({res2['moma_growth']/wt_growth*100:.1f}% retained)")
print("    -> FBA is optimistic: the cell fully reroutes through the PPP;")
print("       MOMA is conservative: fluxes 'resist' large rearrangement,")
print("       growth suffers more - closer to unevolved lab dpgi strains.")

# ============================================================
# 3. 比较敲除株的通量分布差异（以几条核心反应为例）
# ============================================================
print("=" * 60)
print("[3] Flux distribution differences in dpgi: FBA solution vs MOMA solution")
print("    (top 8 reactions by absolute difference):")
fba_flux = res2["_fba_sol"].fluxes
moma_flux = res2["_moma_sol"].fluxes
diff = (fba_flux - moma_flux).abs().sort_values(ascending=False)
show_ids = diff.head(8).index.tolist()
cmp_df = pd.DataFrame({
    "WT_pFBA": wt_pfba.fluxes[show_ids].round(3),
    "dpgi_FBA": fba_flux[show_ids].round(3),
    "dpgi_MOMA": moma_flux[show_ids].round(3),
})
print(cmp_df)
print("  The MOMA solution stays closer to the WT pFBA reference, while the")
print("  FBA solution rearranges fluxes substantially.")

# ============================================================
# 4. 画对比条形图（英文标注）
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))

# 左图：生长率对比
labels = ["Wild type", "dgltA\nFBA", "dgltA\nMOMA", "dpgi\nFBA", "dpgi\nMOMA"]
values = [wt_growth, res1["fba_growth"], res1["moma_growth"],
          res2["fba_growth"], res2["moma_growth"]]
colors = ["#2ca02c", "#d62728", "#ff9896", "#d62728", "#ff9896"]
axes[0].bar(labels, values, color=colors)
axes[0].set_ylabel("Growth rate (1/h)")
axes[0].set_title("Gene knockout growth prediction: FBA vs MOMA")
for i, v in enumerate(values):
    axes[0].text(i, max(v, 0.005), f"{v:.3f}", ha="center", va="bottom", fontsize=9)

# 右图：dpgi 中通量差异最大的反应（FBA解 vs MOMA解）
plot_ids = show_ids[:6]
x = range(len(plot_ids))
width = 0.27
axes[1].bar([i - width for i in x], wt_pfba.fluxes[plot_ids],
            width=width, label="WT pFBA", color="#2ca02c")
axes[1].bar([i for i in x], fba_flux[plot_ids],
            width=width, label="dpgi FBA", color="#d62728")
axes[1].bar([i + width for i in x], moma_flux[plot_ids],
            width=width, label="dpgi MOMA", color="#1f77b4")
axes[1].set_xticks(list(x))
axes[1].set_xticklabels(plot_ids, rotation=30, ha="right")
axes[1].set_ylabel("Flux (mmol/gDW/h)")
axes[1].set_title("dpgi flux distribution: MOMA stays closer to wild type")
axes[1].legend()

plt.tight_layout()
png_path = os.path.join(OUTPUT_DIR, "03_MOMA_knockout_comparison.png")
plt.savefig(png_path, dpi=150)
print("=" * 60)
print("Figure saved:", png_path)