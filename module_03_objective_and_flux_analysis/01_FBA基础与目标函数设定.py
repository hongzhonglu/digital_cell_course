# -*- coding: utf-8 -*-
"""
模块03 脚本01：FBA基础与目标函数设定
================================================
【目的】
1. 理解 FBA（Flux Balance Analysis，通量平衡分析）的核心思想：
   在稳态约束 S·v = 0 与通量上下界 lb <= v <= ub 下，
   通过线性规划（LP）最大化（或最小化）某个目标函数 c·v。
2. 学会查看与更换 COBRA 模型的目标函数（生物量、ATP维持、某产物）。
3. 学会解读 solution 对象：objective_value、fluxes、shadow_prices、reduced_costs。

【用法】python 01_FBA基础与目标函数设定.py
【对应课次】课程5：目标函数与优化原理
"""

# ---------------- 导入依赖库 ----------------
import os
from cobra.io import load_model, load_json_model, read_sbml_model  # 模型加载

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
# e_coli_core 是大肠杆菌核心代谢模型：
# 包含 72 个代谢物、95 个反应、137 个基因，覆盖糖酵解、TCA、戊糖磷酸途径等
model = load_ecoli_core()
BIOMASS = find_biomass_rxn(model).id  # 自动探测生物量反应 ID（两套命名通用）
print("=" * 60)
print("Model:", model)
print("Reactions:", len(model.reactions), " Metabolites:", len(model.metabolites),
      " Genes:", len(model.genes))

# ============================================================
# 1. 查看当前目标函数
# ============================================================
# model.objective 返回一个 sympy 表达式，默认目标为最大化生物量反应（ID 由 find_biomass_rxn 自动探测，
# BiGG 当前版为 BIOMASS_Ecoli_core_w_GAM，textbook 版为 Biomass_Ecoli_core）
print("=" * 60)
print("[1] Default objective (objective expression):")
print(model.objective)
print("Optimization direction:", model.objective.direction)  # 'max' 表示最大化

# ============================================================
# 2. 求解 FBA 并解读 solution 对象
# ============================================================
print("=" * 60)
print("[2] FBA with biomass as the objective:")
solution = model.optimize()  # 等价于求解 LP: max c·v, s.t. S·v=0, lb<=v<=ub
print("Solver status:", solution.status)                     # 'optimal' 表示最优
print("Objective value (max growth rate, 1/h):", round(solution.objective_value, 4))
# fluxes 是一个 pandas Series：每个反应的最优通量（mmol/gDW/h）
print("\nTop 10 reactions by |flux|:")
print(solution.fluxes.abs().sort_values(ascending=False).head(10))
print("\nGlucose uptake EX_glc__D_e:", solution.fluxes["EX_glc__D_e"],
      "(negative = uptake)")
print("O2 uptake EX_o2_e:", solution.fluxes["EX_o2_e"])

# shadow_prices（影子价格）：S·v=0 约束的对偶变量，
# 含义：某代谢物边界增加 1 单位时目标函数的变化率（该代谢物对生长的边际价值）
print("\nTop 5 metabolites by |shadow price| (shadow_prices):")
print(solution.shadow_prices.abs().sort_values(ascending=False).head(5))

# reduced_costs（简约成本）：上下界约束的对偶变量，
# 含义：强制某反应多携带 1 单位通量时目标函数将牺牲多少
print("\nTop 5 reactions by |reduced cost|:")
print(solution.reduced_costs.abs().sort_values(ascending=False).head(5))

# ============================================================
# 3. 更换目标函数：比较不同目标下的优化结果
# ============================================================
print("=" * 60)
print("[3] Switching the objective and comparing objective values:")

# 记录各目标下的优化结果，便于课堂对比
results = {}

# --- 目标 A：最大化生物量（默认，细胞生长）---
model.objective = BIOMASS  # 使用自动探测到的生物量反应 ID
model.objective_direction = "minimize" # 或 "maximize"
sol_a = model.optimize()
results["Biomass (growth)"] = sol_a.objective_value



# --- 目标 B：最大化 ATPM（ATP 维持消耗反应，即"最大化ATP浪费"）---
# ATPM 是一个人工维持反应：ATP + H2O -> ADP + Pi + H，正常解中其值为下界 8.39
# 把目标设为最大化 ATPM，可演示"目标函数决定解"的原理
model.objective = "ATPM"
model.objective_direction = "minimize" # 或 "maximize"
sol_b = model.optimize()
results["ATPM (ATP maintenance)"] = sol_b.objective_value
model.objective.expression


model.objective_direction = "maximize"
sol_b = model.optimize()
results["ATPM (ATP maintenance)"] = sol_b.objective_value


# --- 目标 C：最大化某产物分泌，例如乙酸 EX_ac_e ---
model.objective = "EX_ac_e"
sol_c = model.optimize()
results["Acetate secretion (EX_ac_e)"] = sol_c.objective_value

for name, val in results.items():
    print(f"  Objective = {name:<28s} Optimal value = {val:.4f}")

# 注意：不同目标的量纲不同，数值不可直接互相比较，
# 这里演示的是"换一个目标，LP 就给出另一种最优通量分布"。
# 例如以乙酸为目标时，生物量通量通常降到 0（细胞不生长，全力产酸）
print("\nBiomass flux when optimizing acetate:",
      round(sol_c.fluxes[BIOMASS], 4),
      "(shows the objective directly determines the biological meaning of the solution)")