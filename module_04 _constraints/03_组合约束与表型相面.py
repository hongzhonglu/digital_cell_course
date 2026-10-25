# -*- coding: utf-8 -*-
"""
模块04 脚本03：组合约束与表型相面（PhPP）
================================================
【目的】
1. 理解"表型相面"（Phenotype Phase Plane, PhPP）：
   同时扫描两个环境摄取约束（这里是葡萄糖 × 氧气），
   对每个网格点求解 FBA，得到最大生长率的二维相图。
2. 学会用"production_envelope 思路"手写网格 FBA（不依赖现成函数），
   直观看到约束组合如何划分出不同的"表型相"：
   - 碳限制区（氧充足，生长随葡萄糖线性增加）
   - 氧限制/发酵区（葡萄糖充足但氧不足，被迫发酵）
   - 双底物共同限制过渡区
   - 不可行区（两底物都不足）
3. 讨论：把实验测得的摄取速率作为约束代入模型，落在相面上的具体位置
   即给出唯一预测——"约束越多越准"的原因。

【用法】python 03_组合约束与表型相面.py
【对应课次】课程6：GEMs运行所需的约束
【输出】output/03_glucose_oxygen_phenotype_phase_plane.png
"""

# ---------------- 导入依赖库 ----------------
import os
import warnings
import numpy as np
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

# 网格中葡萄糖=0 时 LP 因 ATPM 维持约束而不可行，cobra 会打印 UserWarning；
# 该情况已在代码中按"不生长"处理，此处静音以保持课堂输出整洁
warnings.filterwarnings("ignore", category=UserWarning, module="cobra")

# ============================================================
# 1. 定义二维扫描网格：葡萄糖摄取 0~10，氧气摄取 0~25
# ============================================================
# 约定：摄取速率取正值（书写方便），赋给反应下界时取负
glc_grid = np.linspace(0, 10, 21)   # 21 个葡萄糖摄取水平
o2_grid = np.linspace(0, 25, 26)    # 26 个氧气摄取水平
print("=" * 60)
print(f"[1] Grid: {len(glc_grid)} glucose levels x {len(o2_grid)} O2 levels"
      f" = {len(glc_grid)*len(o2_grid)} FBA runs")

# ============================================================
# 2. 逐网格点求解 FBA（手写 production_envelope 式扫描）
# ============================================================
# 结果矩阵：行=氧气，列=葡萄糖（便于 imshow/pcolormesh 直接绘图）
growth_map = np.zeros((len(o2_grid), len(glc_grid)))
acetate_map = np.zeros_like(growth_map)  # 同时记录乙酸分泌，辅助划分"发酵相"

with model:  # 整个扫描在 with 块内，结束后自动还原模型
    for i, o2 in enumerate(o2_grid):
        model.reactions.EX_o2_e.lower_bound = -o2
        for j, glc in enumerate(glc_grid):
            model.reactions.EX_glc__D_e.lower_bound = -glc
            sol = model.optimize()
            # LP 不可行时 objective_value 可能为 nan，统一处理
            if sol.status == "optimal":
                growth_map[i, j] = max(sol.objective_value, 0.0)
                acetate_map[i, j] = sol.fluxes["EX_ac_e"]
            else:
                growth_map[i, j] = 0.0
                acetate_map[i, j] = 0.0

print("    Scan finished. A few representative points:")
for o2v, glcv in [(25, 10), (10, 10), (5, 10), (5, 5), (0, 10), (10, 0)]:
    i = np.argmin(abs(o2_grid - o2v))
    j = np.argmin(abs(glc_grid - glcv))
    print(f"    O2={o2_grid[i]:>5.1f}, Glc={glc_grid[j]:>4.1f}"
          f" -> growth = {growth_map[i, j]:.4f} 1/h,"
          f" acetate secretion = {acetate_map[i, j]:.2f}")

# ============================================================
# 3. 简单的"相"划分（教学用规则）
# ============================================================
# 规则：
#   生长率<1e-4        → 不可行/不生长区
#   乙酸分泌>1 且 生长>0 → 氧限制发酵区（溢流代谢）
#   其余               → 呼吸生长区（再按是否接近最大生长率细分碳限制）
phase = np.full(growth_map.shape, "", dtype=object)
phase[growth_map < 1e-4] = "no growth"
phase[(growth_map >= 1e-4) & (acetate_map > 1)] = "O2-limited fermentation"
phase[(growth_map >= 1e-4) & (acetate_map <= 1)] = "respiratory growth"
print("\n[2] Grid-point counts per phenotype phase:")
for name in ["respiratory growth", "O2-limited fermentation", "no growth"]:
    print(f"    {name}: {(phase == name).sum()} points")

# ============================================================
# 4. 绘制表型相面热力图（英文标注）
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(14, 5.5))

# 左图：生长率相面
im0 = axes[0].pcolormesh(glc_grid, o2_grid, growth_map, cmap="viridis",
                         shading="auto")
axes[0].set_xlabel("Glucose uptake rate (mmol/gDW/h)")
axes[0].set_ylabel("O2 uptake rate (mmol/gDW/h)")
axes[0].set_title("Phenotype phase plane: max growth rate (1/h)")
cbar0 = plt.colorbar(im0, ax=axes[0])
cbar0.set_label("Growth rate (1/h)")
# 标注默认条件点（葡萄糖10、氧气≈22 的野生型工作点）
axes[0].plot(10, 21.8, "r*", markersize=18,
             label="default condition (~WT operating point)")
axes[0].legend(loc="lower right")

# 右图：乙酸分泌相面（显示发酵区）
im1 = axes[1].pcolormesh(glc_grid, o2_grid, acetate_map, cmap="OrRd",
                         shading="auto")
axes[1].set_xlabel("Glucose uptake rate (mmol/gDW/h)")
axes[1].set_ylabel("O2 uptake rate (mmol/gDW/h)")
axes[1].set_title("Acetate secretion on the same plane (fermentation phase)")
cbar1 = plt.colorbar(im1, ax=axes[1])
cbar1.set_label("Acetate secretion (mmol/gDW/h)")

plt.tight_layout()
png_path = os.path.join(OUTPUT_DIR, "03_glucose_oxygen_phenotype_phase_plane.png")
plt.savefig(png_path, dpi=150)
print("\nFigure saved:", png_path)

# ============================================================
# 5. 讨论打印
# ============================================================
print("=" * 60)
print("[3] Discussion: how constraints sharpen predictions")
print("  1) Saying only 'glucose is present' lets the model pick ANY point on")
print("     the plane - the prediction is a whole region;")
print("  2) Plugging measured uptake rates in as constraints pins down ONE")
print("     operating point - the more constraints, the sharper the prediction;")
print("  3) The plane also explains experiments: high glucose + low O2 ->")
print("     acetate overflow (fermentation phase), the metabolic state that")
print("     industrial E. coli fermentations try to avoid.")
print("Script 03 finished.")
