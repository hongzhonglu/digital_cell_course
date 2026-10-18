# -*- coding: utf-8 -*-
# =====================================================================
# 01_从零构建toy模型.py
# ---------------------------------------------------------------------
# 模块二《GEMs构建与表型预测》（对应课程4）第 1 个上机脚本。
#
# 目的：不用现成模型，用 cobra.Model / Reaction / Metabolite 手工构建一个
#       “最小糖酵解” toy 模型，理解 GEM 的三要素（代谢物、反应、边界条件），
#       再添加交换反应与生物量反应，设定目标函数并用 FBA 求解。
#
# Toy 模型途径（简化版糖酵解 + 呼吸 + 生长）：
#     葡萄糖(glc) --v1--> 6-磷酸果糖(f6p) --v2--> 3-磷酸甘油醛(g3p)
#     g3p --v3(产能)--> 丙酮酸(pyr) ，pyr --v4--> 排出
#     代谢物前体(f6p+g3p) --biomass--> 生物量（目标函数）
#
# 用法：python 01_从零构建toy模型.py
# 依赖：cobra
# =====================================================================
"""Build a minimal glycolysis toy model from scratch with cobra and solve it by FBA."""

from cobra import Metabolite, Model, Reaction

# ----------------------------------------------------------------------
# 1. 创建空模型与代谢物 Metabolite
# ----------------------------------------------------------------------
toy = Model("toy_glycolysis")          # 新建一个空模型，id 为 toy_glycolysis

# Metabolite 的关键属性：id、化学式(formula)、名称(name)、所在区室
# （c=细胞质，e=细胞外）；formula 可留空，但写上更规范
glc_e = Metabolite("glc_e", formula="C6H12O6", name="Glucose (extracellular)", compartment="e")
glc_c = Metabolite("glc_c", formula="C6H12O6", name="Glucose", compartment="c")
f6p_c = Metabolite("f6p_c", formula="C6H13O9P", name="Fructose-6-phosphate", compartment="c")
g3p_c = Metabolite("g3p_c", formula="C3H7O6P", name="Glyceraldehyde-3-phosphate", compartment="c")
pyr_c = Metabolite("pyr_c", formula="C3H3O3", name="Pyruvate", compartment="c")
pyr_e = Metabolite("pyr_e", formula="C3H3O3", name="Pyruvate (extracellular)", compartment="e")
atp_c = Metabolite("atp_c", formula="C10H12N5O13P3", name="ATP", compartment="c")
adp_c = Metabolite("adp_c", formula="C10H12N5O10P2", name="ADP", compartment="c")

# ----------------------------------------------------------------------
# 2. 定义一个辅助函数，批量创建反应，减少重复代码
# ----------------------------------------------------------------------
def make_reaction(rxn_id, name, stoichiometry, lb=0.0, ub=1000.0):
    """Create a cobra.Reaction and add metabolites with stoichiometric coefficients.

    Parameters
    ----------
    rxn_id : str
        Reaction identifier.
    name : str
        Human-readable reaction name.
    stoichiometry : dict
        {Metabolite: coefficient}; negative = substrate, positive = product.
    lb, ub : float
        Lower/upper flux bounds (mmol/gDW/h); lb < 0 makes the reaction reversible.

    Returns
    -------
    cobra.Reaction
        The constructed reaction object.
    """
    reaction = Reaction(rxn_id)             # 新建反应
    reaction.name = name                    # 设置名称
    reaction.lower_bound = lb               # 通量下界
    reaction.upper_bound = ub               # 通量上界
    reaction.add_metabolites(stoichiometry) # 一次性添加底物与产物
    return reaction


# ----------------------------------------------------------------------
# 3. 逐步搭建途径：运输 -> 糖酵解 -> 排出
# ----------------------------------------------------------------------
# v0：葡萄糖从胞外运输到胞内（glc_e -> glc_c）
v0 = make_reaction("GLCt", "Glucose transport",
                   {glc_e: -1.0, glc_c: 1.0})

# v1：葡萄糖磷酸化并异构为 6-磷酸果糖（消耗 1 ATP），对应 HXK+PGI
v1 = make_reaction("v1", "Glucose -> fructose-6-phosphate",
                   {glc_c: -1.0, atp_c: -1.0, f6p_c: 1.0, adp_c: 1.0})

# v2：f6p 裂解为 2 分子 g3p（合并 PFK+FBA 两步）
v2 = make_reaction("v2", "Fructose-6-phosphate -> 2 glyceraldehyde-3-phosphate",
                   {f6p_c: -1.0, g3p_c: 2.0})

# v3：g3p 氧化为丙酮酸并产生 2 ATP（合并下游糖酵解，净产能）
v3 = make_reaction("v3", "Glyceraldehyde-3-phosphate -> pyruvate + 2 ATP",
                   {g3p_c: -1.0, adp_c: -2.0, pyr_c: 1.0, atp_c: 2.0})

# v4：丙酮酸输出到胞外（pyr_c -> pyr_e）
v4 = make_reaction("PYRt", "Pyruvate transport",
                   {pyr_c: -1.0, pyr_e: 1.0})

# ----------------------------------------------------------------------
# 4. 添加生物量反应（目标函数）与 ATP 维持消耗
# ----------------------------------------------------------------------
# 生物量反应：消耗 1 f6p + 1 g3p + 4 ATP 产生 1 单位“生物量”
biomass = make_reaction("BIOMASS", "Biomass synthesis (objective)",
                        {f6p_c: -1.0, g3p_c: -1.0,
                         atp_c: -4.0, adp_c: 4.0})

# ATP 维持反应：atp -> adp，模拟细胞非生长相关的能量消耗
atpm = make_reaction("ATPM", "ATP maintenance",
                     {atp_c: -1.0, adp_c: 1.0}, lb=1.0)  # lb=1 强制至少消耗 1 单位 ATP

# ----------------------------------------------------------------------
# 5. 添加交换反应（boundary），定义培养基：允许摄入葡萄糖、排出丙酮酸
# ----------------------------------------------------------------------
ex_glc = make_reaction("EX_glc_e", "Glucose exchange",
                       {glc_e: -1.0}, lb=-10.0)   # 下界 -10：最多摄入 10
ex_pyr = make_reaction("EX_pyr_e", "Pyruvate exchange",
                       {pyr_e: -1.0}, lb=0.0)     # 只允许排出，不允许摄入

# 把所有反应加入模型
toy.add_reactions([v0, v1, v2, v3, v4, biomass, atpm, ex_glc, ex_pyr])

# 给反应补充 GPR 与 subsystem 信息（真实模型中用于基因敲除模拟与途径分类）
# 写法示例："a and b" 表示需要两个基因共同编码的酶复合体；"a or b" 表示同工酶
v1.gene_reaction_rule = "glk"           # 单个基因催化（虚构基因名，仅演示写法）
v2.gene_reaction_rule = "(pfkA or pfkB)"  # 同工酶逻辑
v3.gene_reaction_rule = "gapA and pgk"    # 酶复合体/多步逻辑
for rxn in [v0, v1, v2, v3, v4]:
    rxn.subsystem = "Glycolysis"          # subsystem 标注反应所属途径
print("GPR of reaction v2:", v2.gene_reaction_rule,
      "| parsed genes:", sorted(g.id for g in v2.genes))

print("\n===== Toy model built =====")
print("Model:", toy)
print("Reactions:", len(toy.reactions), "| Metabolites:", len(toy.metabolites),
      "| Genes:", len(toy.genes))
print("\nAll reaction equations and bounds:")
for rxn in toy.reactions:
    print(f"  {rxn.id:10s} [{rxn.lower_bound:>6}, {rxn.upper_bound:>6}]  {rxn.reaction}")

# model.medium 汇总当前培养基（下界 < 0 的交换反应的摄取上限）
print("\nmodel.medium:", toy.medium)

# ----------------------------------------------------------------------
# 6. 设定目标函数并求解 FBA
# ----------------------------------------------------------------------
toy.objective = "BIOMASS"            # 目标：最大化生物量反应通量
print("\nObjective function:", toy.objective.expression)

solution = toy.optimize()            # 求解线性规划
print("\n===== FBA result =====")
print("Solver status:", solution.status)
print(f"Max growth rate (biomass flux): {solution.objective_value:.4f}")
print("\nOptimal flux of each reaction:")
for rxn_id, flux in solution.fluxes.items():
    print(f"  {rxn_id:10s} flux = {flux:8.3f}")

# model.summary() 给出摄取/分泌/目标的一览表，适合课堂快速检查
print("\n===== model.summary() =====")
toy.summary()

# ----------------------------------------------------------------------
# 7. 课堂讨论：把葡萄糖摄取上限改成 5，生长速率如何变化？
# ----------------------------------------------------------------------
toy.reactions.EX_glc_e.lower_bound = -5.0    # 限制葡萄糖供应减半
sol2 = toy.optimize()
print(f"\nAfter limiting glucose uptake to 5: max growth = {sol2.objective_value:.4f} "
      f"(was {solution.objective_value:.4f})")

print("\n=== Script 01_从零构建toy模型.py finished ===")
