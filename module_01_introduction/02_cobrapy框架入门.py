# -*- coding: utf-8 -*-
# =====================================================================
# 02_cobrapy框架入门.py
# ---------------------------------------------------------------------
# 模块一《入门知识介绍》（对应课程3）第 2 个上机脚本。
#
# 目的：认识 COBRApy 的核心对象 Model / Reaction / Metabolite / Gene，
#       遍历大肠杆菌核心模型 e_coli_core，查看反应式、上下界、GPR 规则，
#       并运行一次 FBA（Flux Balance Analysis，通量平衡分析）解读结果。
#
# 用法：python 02_cobrapy框架入门.py
# 依赖：cobra（pip install cobra）
# =====================================================================
"""Guided tour of the COBRApy core objects using the e_coli_core model."""

# model_loader 会优先在线加载 BiGG 模型，无网络时自动使用本地缓存
# 在左侧的项目目录中，找到你的文件夹 module_01_introduction。   右键点击 该文件夹 module_01_introduction。   在弹出的菜单中，找到并点击 "Mark Directory as" -> "Sources Root"（将其标记为蓝色）。
from model_loader import load_ecoli_core

# ----------------------------------------------------------------------
# 1. 加载模型并查看总体规模
# ----------------------------------------------------------------------
model = load_ecoli_core()          # 加载 e_coli_core 模型（cobra.Model 对象）
print("\n===== Basic model information =====")
print("Model object:", model)      # print 模型会显示摘要
print("Model ID:", model.id)
print("Number of reactions:", len(model.reactions))      # Reaction 个数
print("Number of metabolites:", len(model.metabolites))  # Metabolite 个数
print("Number of genes:", len(model.genes))              # Gene 个数

# ----------------------------------------------------------------------
# 2. 遍历反应 Reaction：查看反应式、上下界
# ----------------------------------------------------------------------
print("\n===== Details of the first 5 reactions =====")
for reaction in model.reactions[:5]:           # DictList 支持切片，取前 5 个
    # reaction.reaction 属性给出化学计量反应式（含系数与方向箭头）
    # reaction.lower_bound / upper_bound 为通量下界/上界（mmol/gDW/h）
    print(f"Reaction {reaction.id:12s} | {reaction.reaction}")
    print(f"{'':15s} lower_bound={reaction.lower_bound:>8}, upper_bound={reaction.upper_bound:>8}")

# 通过 id 快速取出一个反应（以 ATP 合酶 ATPS4r 为例）
atps = model.reactions.get_by_id("ATPS4r")
print("\nATPS4r equation:", atps.reaction)
print("ATPS4r name:", atps.name)               # name 为可读的反应名称
print("ATPS4r reversibility:", "reversible" if atps.reversibility else "irreversible")
# bounds 属性返回 (下界, 上界) 元组
print("ATPS4r flux bounds:", atps.bounds)
# annotation 属性给出该反应在其它数据库中的对应编号（如 BiGG、KEGG、EC 号）
print("ATPS4r annotation keys:", list(atps.annotation.keys()))

# 查看反应参与的代谢物及其化学计量系数（字典：Metabolite -> 系数）
print("ATPS4r metabolites and coefficients:")
for met, coef in atps.metabolites.items():
    print(f"   {met.id:10s} coefficient = {coef}")   # 负数为底物，正数为产物

# query() 可按名称模糊搜索反应，例如找出所有名字含 "biomass" 的反应
print("\nReactions whose name contains 'biomass':")
for rxn in model.reactions.query("biomass", attribute="name"):
    print(f"   {rxn.id:30s} | {rxn.name}")

# ----------------------------------------------------------------------
# 3. 遍历代谢物 Metabolite：查看名称、化学式、所属区室
# ----------------------------------------------------------------------
print("\n===== First 5 metabolites =====")
for met in model.metabolites[:5]:
    # met.compartment 为区室：c=细胞质, e=细胞外, p=周质
    print(f"{met.id:12s} | name: {met.name:25s} | formula: {met.formula} | compartment: {met.compartment}")

# 统计某代谢物参与的反应数量（以胞内 ATP 为例）
atp = model.metabolites.get_by_id("atp_c")
print(f"\nNumber of reactions involving atp_c: {len(atp.reactions)}")

# ----------------------------------------------------------------------
# 4. 基因 Gene 与 GPR 规则（Gene-Protein-Reaction association）
# ----------------------------------------------------------------------
print("\n===== GPR (gene-protein-reaction) rules =====")
# GPR 用 and/or 逻辑描述“哪些基因编码的酶催化该反应”
for rxn_id in ["PFK", "GAPD", "ATPS4r"]:
    rxn = model.reactions.get_by_id(rxn_id)
    print(f"Reaction {rxn_id:8s} GPR: {rxn.gene_reaction_rule}")
    # reaction.genes 自动解析出涉及的基因对象集合
    print(f"{'':11s}associated genes: {sorted(g.id for g in rxn.genes)}")

# 反过来，查看某个基因参与了哪些反应
gene = model.genes.get_by_id("b1779")          # gapA，编码 GAPD 的酶
print(f"\nReactions associated with gene b1779: {[r.id for r in gene.reactions]}")

# 演示单基因敲除：gene.knock_out() 会关闭该基因 GPR 控制的所有反应
# with model: 上下文保证退出后模型自动还原，不影响后续分析
wt_growth = model.optimize().objective_value
with model:
    model.genes.get_by_id("b1779").knock_out()     # 敲除 gapA
    ko_growth = model.optimize().objective_value   # 敲除后再求解
print(f"Gene knock-out demo (b1779/gapA): WT growth = {wt_growth:.4f}, "
      f"knock-out growth = {ko_growth:.4f} 1/h")

# ----------------------------------------------------------------------
# 5. 查看边界反应（exchange）——模型的“培养基”设置
# ----------------------------------------------------------------------
print("\n===== Exchange reactions (medium composition) =====")
# model.exchanges 给出所有交换反应；下界 < 0 表示允许摄入
for ex in model.exchanges:
    if ex.lower_bound < 0:                     # 只打印允许摄入的营养物
        print(f"{ex.id:15s} lower_bound={ex.lower_bound:>8}  equation: {ex.reaction}")

# model.medium 直接以字典形式给出当前培养基中各碳源/营养物的摄取上限
print("\nmodel.medium (uptake limits):", model.medium)

# ----------------------------------------------------------------------
# 6. 运行一次 FBA 并解读结果
# ----------------------------------------------------------------------
print("\n===== Running FBA =====")
print("Objective function:", model.objective.expression)  # 默认目标为最大化生物量反应
solution = model.optimize()                    # 求解线性规划，返回 Solution 对象
print("Solver status:", solution.status)       # "optimal" 表示求得最优解
# objective_value 即最大生长速率（单位 1/h）
print(f"Max growth rate: {solution.objective_value:.4f} 1/h")

# solution.fluxes 是 pandas.Series，索引为反应 id，值为最优通量
print("\nKey reaction fluxes (mmol/gDW/h):")
for rxn_id in ["EX_glc__D_e",   # 葡萄糖摄取
               "EX_o2_e",       # 氧气摄取
               "EX_co2_e",      # 二氧化碳排出
               "EX_ac_e",       # 乙酸排出
               "BIOMASS_Ecoli_core_w_GAM"]:  # 生物量（生长）
    flux = solution.fluxes[rxn_id]
    print(f"  {rxn_id:28s} flux = {flux:8.3f}")

# 取通量绝对值最大的 10 个反应，展示细胞内“最繁忙”的途径
top_flux = solution.fluxes.abs().sort_values(ascending=False).head(10)
print("\nTop 10 reactions by |flux|:")
for rxn_id in top_flux.index:
    print(f"  {rxn_id:12s} flux = {solution.fluxes[rxn_id]:8.3f}")

# solution.to_frame() 把结果转成 DataFrame，方便筛选（如只看有通量的反应）
flux_frame = solution.to_frame()               # 列为 fluxes / reduced_costs
active = flux_frame[flux_frame["fluxes"].abs() > 1e-4]
print(f"\nReactions with |flux| > 1e-4: {len(active)} / {len(flux_frame)} "
      "(use solution.to_frame() for further pandas analysis)")

# model.summary() 一键打印模型的摄取/分泌/目标概览（课堂快速检查用）
print("\n===== model.summary() overview =====")
model.summary()

print("\n=== Script 02_cobrapy框架入门.py finished ===")
