# -*- coding: utf-8 -*-
"""
模块03 脚本05：通量取样（Flux Sampling）——均匀刻画整个解空间
================================================
【目的】
1. 掌握 cobra.sampling（ACHRSampler / OptGPSampler）的用法：
   在稳态约束 S·v=0 与通量边界定义的"可行解多面体"内做均匀随机取样，
   得到一批（如 1000 个）可行通量分布，用统计方法刻画解空间。
2. 理解"取样"与 FVA 的互补关系：
   - FVA 给出每个反应"独立"的通量区间（区间的边界），但无法回答
     "两个反应能不能同时取到各自区间端点"——即反应间的相关性；
   - 取样给出的是一批"完整"的通量向量，可以直接计算任意反应间的
     相关/互斥关系（散点图、相关系数）以及每个反应的边际分布（直方图）。
   一句话：FVA 看"边缘"，取样看"形状"。
3. 实操：用 OptGPSampler 取样 1000 个点，分析
   Biomass / EX_ac_e / EX_co2_e 的通量分布直方图与两两相关性；
   结果存 png（直方图+散点矩阵）与 csv（样本表 + 相关系数表）。

【用法】python 05_flux_sampling.py（整体运行约 0.5~2 分钟，取样最耗时）
【对应课次】课程5：目标函数与优化原理
【输出】
  output/05_flux_sampling_distributions.png  直方图 + 相关性散点图
  output/05_flux_samples.csv                 1000 个通量样本（关键反应）
  output/05_flux_correlation.csv             关键反应 Pearson 相关系数表
"""

# ---------------- 导入依赖库 ----------------
import os
import time
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # 无界面环境下保存图片
import matplotlib.pyplot as plt
from cobra.io import load_model, load_json_model, read_sbml_model  # 模型加载
from cobra.sampling import OptGPSampler, ACHRSampler               # 两种取样器

# 图输出目录
OUTPUT_DIR = os.path.join(os.path.dirname(os.path.abspath(__file__)), "output")
os.makedirs(OUTPUT_DIR, exist_ok=True)

# ---------------- 取样参数（在运行时间与样本质量之间折中） ----------------
N_SAMPLES = 1000   # 取样点数：1000 足够看清分布形状与相关性
THINNING = 100     # 链上每 100 步保留 1 个样本（减小相邻样本的自相关）
                   # OptGP 的 warmup 点数默认为 2*len(reactions)，足够小模型使用；
                   # 大模型可给 OptGPSampler(..., thinning=...) 并减小 batch 以省时


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


# ---------------- 加载模型 ----------------
model = load_ecoli_core()
BIOMASS = find_biomass_rxn(model).id  # 自动探测生物量反应 ID（两套命名通用）
opt_growth = model.optimize().objective_value
print("=" * 60)
print(f"Model loaded. Optimal growth rate (for reference) = {opt_growth:.4f} 1/h")
print("Note: sampling is done on the UNCONSTRAINED objective (no minimum")
print("growth enforced) - samples cover the whole feasible steady-state space.")

# ============================================================
# 1. 构造 OptGPSampler 并取样
# ============================================================
# OptGP（Optimal Gaussian Perturbation）取样器：
#   - 先自动生成一批 warmup 点（默认 2*反应数 个，分散在可行多面体中），
#   - 再沿"扰动的方向"走 hit-and-run 链，理论上能更均匀地覆盖大空间。
# ACHRSampler（Artificial Centering Hit-and-Run）是另一种经典实现：
#   - 无需 LP 求 warmup 的精细度要求低，但在狭长空间中混合较慢；
#   两者用法完全一致：sampler = ACHRSampler(model, thinning=THINNING)。
# processes=1：课堂演示用单进程（多进程 pickle 模型反而可能更慢）
print("=" * 60)
print(f"[1] Building OptGPSampler (thinning={THINNING}, processes=1) ...")
t0 = time.time()
sampler = OptGPSampler(model, processes=1, thinning=THINNING)
print(f"    Sampler ready in {time.time()-t0:.1f} s"
      f" (including {sampler.n_warmup} warmup points).")

# 取样 N_SAMPLES 个可行通量分布：返回 N x 95 的 DataFrame，列为反应 ID
t0 = time.time()
samples = sampler.sample(N_SAMPLES)
print(f"    Drew {samples.shape[0]} samples x {samples.shape[1]} reactions"
      f" in {time.time()-t0:.1f} s.")

# 校验：validate 用 LP 逐个检查样本是否严格可行（'v' = valid）
# 教学上演示其用法；大批量校验较慢，这里只抽查前 50 个
check = pd.Series(sampler.validate(samples.iloc[:50]))
print("    Validity check on the first 50 samples:",
      check.value_counts().to_dict(), "('v' = valid)")

# ============================================================
# 2. 关键反应的通量分布（边际分布）
# ============================================================
key_ids = [BIOMASS, "EX_glc__D_e", "EX_o2_e", "EX_ac_e", "EX_co2_e"]
key_ids = [r for r in key_ids if r in samples.columns]  # 防御：只保留存在的反应
key_samples = samples[key_ids]

print("=" * 60)
print("[2] Marginal distributions of key reactions (from samples):")
print(key_samples.describe().round(3))

# 对照：FBA 最优解给出的"单点"通量，与取样的均值/区间对比
opt_flux = model.optimize().fluxes[key_ids]
cmp = pd.DataFrame({
    "FBA optimum": opt_flux.round(3),
    "sample mean": key_samples.mean().round(3),
    "sample min": key_samples.min().round(3),
    "sample max": key_samples.max().round(3),
})
print("\n    FBA single-point solution vs sampled distribution:")
print(cmp)
print("    -> The FBA optimum sits at the EDGE of the sampled range;")
print("       the 'typical' (mean) flux can be far from the optimum.")

# ============================================================
# 3. 反应间相关性：取样独有的信息（FVA 给不出）
# ============================================================
corr = key_samples.corr()  # Pearson 相关系数矩阵
print("=" * 60)
print("[3] Pearson correlation between key reactions (samples only):")
print(corr.round(3))
print("    Reading: strongly negative/positive pairs reveal trade-offs or")
print("    couplings (e.g. acetate secretion vs growth = overflow trade-off).")
print("    FVA CANNOT give this: it varies each reaction independently.")

# ============================================================
# 4. 保存结果（csv）
# ============================================================
csv_samples = os.path.join(OUTPUT_DIR, "05_flux_samples.csv")
csv_corr = os.path.join(OUTPUT_DIR, "05_flux_correlation.csv")
key_samples.round(6).to_csv(csv_samples, index=False)
corr.round(4).to_csv(csv_corr, index_label="reaction")
print("=" * 60)
print("Samples saved:", csv_samples)
print("Correlation table saved:", csv_corr)

# ============================================================
# 5. 绘图：直方图（边际分布）+ 散点矩阵（相关性）（英文标注）
# ============================================================
# 布局：len(key_ids) x len(key_ids) 网格——对角线为直方图，非对角线为散点图
n = len(key_ids)
fig, axes = plt.subplots(n, n, figsize=(2.6 * n, 2.4 * n))
for i, ri in enumerate(key_ids):
    for j, rj in enumerate(key_ids):
        ax = axes[i, j]
        if i == j:
            # 对角线：该反应的通量分布直方图
            ax.hist(key_samples[ri], bins=30, color="#1f77b4", alpha=0.85)
            # 叠加 FBA 最优解位置（红竖线），直观看到"最优≠典型"
            ax.axvline(opt_flux[ri], color="red", ls="--", lw=1.5,
                       label="FBA optimum" if i == 0 else None)
        else:
            # 非对角线：两反应的通量散点（揭示相关性）
            ax.scatter(key_samples[rj], key_samples[ri],
                       s=4, alpha=0.3, color="#2ca02c")
            ax.text(0.05, 0.9, f"r = {corr.loc[ri, rj]:.2f}",
                    transform=ax.transAxes, fontsize=9, color="black")
        if i == n - 1:
            ax.set_xlabel(rj, fontsize=8)
        if j == 0:
            ax.set_ylabel(ri, fontsize=8)
        ax.tick_params(labelsize=7)

fig.suptitle(f"Flux sampling ({N_SAMPLES} samples): histograms (diagonal)"
             " and pairwise scatter (off-diagonal)", y=1.0)
axes[0, 0].legend(fontsize=8)
plt.tight_layout()
png_path = os.path.join(OUTPUT_DIR, "05_flux_sampling_distributions.png")
plt.savefig(png_path, dpi=150)
print("Figure saved:", png_path)

# ============================================================
# 6. 小结
# ============================================================
print("=" * 60)
print("[Summary] Flux sampling vs FVA: complementary views of the solution space")
print("  FVA      -> independent min/max range per reaction (the 'edges');")
print("  Sampling -> full flux vectors drawn uniformly from the space:")
print("              marginal distributions (histograms), joint behaviour")
print("              (scatter, correlations), and 'typical' vs 'optimal' fluxes.")
print("  Use FVA for quick range/confidence checks; use sampling when the")
print("  JOINT distribution of fluxes (trade-offs, couplings) matters.")
print("Script 05 finished.")
