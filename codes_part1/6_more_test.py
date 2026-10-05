import cobra
import matplotlib.pyplot as plt


# module 1 模型加载与基本统计
# 加载模型（假设 iJO1366.xml 在当前目录）
model = cobra.io.read_sbml_model("data/iJO1366.xml")

# 基本统计
print(f"模型名称: {model.id}")
print(f"反应总数: {len(model.reactions)}")
print(f"代谢物总数: {len(model.metabolites)}")
print(f"基因总数: {len(model.genes)}")
print(f"生物质反应 (Biomass): {model.objective}")

# 展示前5个反应
print("\n前5个反应:")
for rxn in list(model.reactions)[:5]:
    print(f"- {rxn.id}: {rxn.name}")

# 可视化：反应与代谢物数量分布（简单柱状图）
fig, ax = plt.subplots(1, 1, figsize=(8, 5))
categories = ['Reactions', 'Metabolites', 'Genes']
counts = [len(model.reactions), len(model.metabolites), len(model.genes)]
ax.bar(categories, counts, color=['blue', 'green', 'red'])
ax.set_ylabel('Count')
ax.set_title('iJO1366 Model Overview')
plt.show()



# Module 2 FBA分析
# 运行 FBA（默认最小介质：葡萄糖 + 氧气）
solution = model.optimize()

# 输出生长率（生物质通量）
print(f"最大生长率 (μ): {solution.objective_value:.4f} h⁻¹")

# 展示高通量反应
print("\n高通量反应 (前5):")
flux_df = cobra.flux_analysis.flux_variability_analysis(model).T
high_flux = solution.fluxes.abs() > 10
for rxn_id, flux in solution.fluxes[high_flux].iteritems():
    print(f"- {rxn_id}: {flux:.2f}")

# 可视化：通量分布直方图
plt.figure(figsize=(8, 5))
plt.hist(solution.fluxes, bins=50, alpha=0.7, edgecolor='black')
plt.xlabel('Flux Value')
plt.ylabel('Frequency')
plt.title('Flux Distribution in iJO1366 (FBA)')
plt.axvline(0, color='red', linestyle='--', label='Zero Flux')
plt.legend()
plt.show()


# module 3 单基因敲除分析
from cobra.flux_analysis import single_gene_deletion

# 单基因敲除分析（前50个基因，随机选择或指定）
genes_to_test = [gene.id for gene in list(model.genes)[:200]]  # 或指定如 'b0001'


deletion_results = single_gene_deletion(model, gene_list=genes_to_test)

# 输出：敲除后生长率变化
print("基因敲除结果:")
print(deletion_results)

# 可视化：生长率保留比例
growth_retention = deletion_results['growth'] / solution.objective_value
plt.figure(figsize=(10, 6))
plt.bar(range(len(genes_to_test)), growth_retention, color='orange', alpha=0.7)
plt.xlabel('Gene Index')
plt.ylabel('Growth Retention Ratio')
plt.title('Single Gene Deletion Effects on Growth')
plt.xticks(range(len(genes_to_test)), genes_to_test, rotation=45)
plt.axhline(1, color='red', linestyle='--', label='Wild Type')
plt.legend()
plt.tight_layout()
plt.show()

# 检查必需基因（生长率为0）
essential = deletion_results[deletion_results['growth'] <= 0.01]
if not essential.empty:
    print(f"\n潜在必需基因: {essential.index.tolist()}")




# module 4 可变性分析
from cobra.flux_analysis import flux_variability_analysis

# 优化最小介质（使用 OptKnock-like 简单实现：最小化必需营养）
# 步骤1: 关闭所有交换反应
#for rxn in model.exchanges:
#    rxn.bounds = (0, 0)  # 无摄取/分泌
model = cobra.io.read_sbml_model("data/iJO1366.xml")

# 步骤2: 打开必需交换（葡萄糖、氧、NH4+ 等）
model.reactions.EX_glc__D_e.bounds = (-10, 0)  # 葡萄糖摄取
model.reactions.EX_o2_e.bounds = (-20, 0)     # 氧气
model.reactions.EX_nh4_e.bounds = (-20, 0)    # 铵离子

# 运行 FBA
solution_opt = model.optimize()
print(f"优化介质下生长率: {solution_opt.objective_value:.4f} h⁻¹")

# FVA：分析关键反应通量变化（e.g., 糖酵解通路）
target_reactions = [rxn.id for rxn in model.reactions if 'PGI' in rxn.id or 'PFK' in rxn.id]

allrxn = []
for rxn in model.reactions:
    print(rxn.id)
    allrxn.append(rxn.id)

fva_results = flux_variability_analysis(model, reaction_list=target_reactions)

print("\nFVA 结果 (min-max flux):")
print(fva_results)

# 可视化：FVA 箱线图
fva_results.plot(kind='bar', y=['minimum', 'maximum'], figsize=(10, 6))
plt.title('Flux Variability in Key Reactions')
plt.ylabel('Flux Bounds')
plt.xticks(rotation=45)
plt.show()


# module 5 预测底物使用
import cobra
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# 加载模型
model = cobra.io.read_sbml_model("data/iJO1366.xml")

# 定义碳源列表（交换反应 ID，从 BiGG 数据库查询）
carbon_sources = {
    'Glucose': 'EX_glc__D_e',
    'Succinate': 'EX_succ_e',
    'Lactate': 'EX_lac__D_e',
    'Acetate': 'EX_ac_e'
}

# 通用设置：最小介质（氧、氮等固定）
#    if rxn.id not in ['EX_o2_e', 'EX_nh4_e', 'EX_pi_e', 'EX_so4_e', 'EX_h2o_e']:
#        rxn.bounds = (0, 0)  # 关闭非必需交换
model.reactions.EX_o2_e.bounds = (-20, 0)  # 需氧
model.reactions.EX_nh4_e.bounds = (-20, 0)
model.reactions.EX_pi_e.bounds = (-20, 0)
model.reactions.EX_so4_e.bounds = (-20, 0)
model.reactions.EX_h2o_e.bounds = (-1000, 1000)  # 水自由

# 预测函数：为每个碳源设置摄取率，运行 FBA
results = []
with model:
    for name, rxn_id in carbon_sources.items():
        model.reactions.get_by_id(rxn_id).bounds = (-10, 0)  # 固定摄取 10 mmol/gDW/h
        solution = model.optimize()
        growth = solution.objective_value
        yield_biomass = growth / abs(model.reactions.get_by_id(rxn_id).flux) if growth > 0 else 0
        results.append({'Substrate': name, 'Growth Rate': growth, 'Biomass Yield': yield_biomass})
        model.reactions.get_by_id(rxn_id).bounds = (0, 0)  # 重置

# 输出结果
df_results = pd.DataFrame(results)
print("底物使用预测结果:")
print(df_results.round(4))

# 可视化：生长率柱状图
plt.figure(figsize=(8, 5))
plt.bar(df_results['Substrate'], df_results['Growth Rate'], color='skyblue', alpha=0.7)
plt.ylabel('Growth Rate (h⁻¹)')
plt.title('Predicted Growth on Different Substrates')
plt.xticks(rotation=45)
plt.show()



# module 6 相平面分析

from cobra.flux_analysis import production_envelope

# 加载模型（继承模块1设置，或重置）
model = cobra.io.read_sbml_model("data/iJO1366.xml")
# 基本介质设置（同模块1）
BIOMASS_ID = 'BIOMASS_Ec_iJO1366_core_53p95M'
# 定义产物：琥珀酸交换反应
product_rxn = model.reactions.EX_succ_e
product_rxn.bounds = (0, 1000)  # 允许分泌
growth_rxn = model.reactions.get_by_id(BIOMASS_ID)

model.reactions.EX_glc__D_e.bounds = (-10, 0)  # 葡萄糖摄取


# 步骤3: 设置 objective 为产物（最大化琥珀酸），变异生长通量
original_objective = model.objective  # 保存原目标
model.objective = product_rxn  # 临时设为产物

envelope = production_envelope(
    model,
    reactions=[growth_rxn],  # 只变异生长反应（从 0 到 max）
    points=20  # 采样点
)

# 重置 objective
model.objective = original_objective


# 输出关键点
print("相平面关键点 (Growth vs. Succinate Yield):")
print(envelope.round(4))

# 可视化：Pareto 前沿
plt.figure(figsize=(8, 6))
plt.plot(envelope['BIOMASS_Ec_iJO1366_core_53p95M'], envelope['flux_maximum'], 'o-', color='red')
plt.xlabel('Growth Rate (h⁻¹)')
plt.ylabel('Succinate Flux (mmol/gDW/h)')
plt.title('Phase Plane: Growth vs. Succinate Production')
plt.grid(True, alpha=0.3)
# 添加最大产量点
max_yield_point = envelope.loc[envelope['flux_maximum'].idxmax()]
plt.annotate(f'Max Yield: {max_yield_point["flux_maximum"]:.2f}',
             (max_yield_point['BIOMASS_Ec_iJO1366_core_53p95M'], max_yield_point['flux_maximum']),
             arrowprops=dict(arrowstyle='->'))
plt.show()



# module 7 通量采样分析
from cobra.sampling import sample

# 加载模型（继承介质设置）
model = cobra.io.read_sbml_model("data/iJO1366.xml")
# 基本介质（同模块1）

# 运行通量取样：Hit-and-run 采样，生成 1000 个样本
samples = sample(model, n=1000, seed=42)  # n=样本数，seed=可重复

# 输出：样本统计（e.g., 关键反应均值/方差）
key_rxns = ['PGI', 'PFK']  # 糖酵解关键酶
print("通量取样统计 (均值 ± 标准差):")
for rxn_id in key_rxns:
    if rxn_id in samples.columns:
        mean_flux = samples[rxn_id].mean()
        std_flux = samples[rxn_id].std()
        print(f"- {rxn_id}: {mean_flux:.4f} ± {std_flux:.4f}")

# 可视化：采样分布直方图（示例一个反应）
rxn_sample = 'PGI'
plt.figure(figsize=(10, 4))
plt.subplot(1, 2, 1)
plt.hist(samples[rxn_sample], bins=50, alpha=0.7, edgecolor='black')
plt.xlabel(f'{rxn_sample} Flux')
plt.ylabel('Frequency')
plt.title('Flux Sampling Distribution')

# 子图2：相关性热图（前5反应）
top_corr = samples[key_rxns[:2]].corr()
plt.subplot(1, 2, 2)
import seaborn as sns
sns.heatmap(top_corr, annot=True, cmap='coolwarm', center=0)
plt.title('Flux Correlations')
plt.tight_layout()
plt.show()





# module 8 自定义目标函数
from cobra import Model, Reaction

# 加载模型
model = cobra.io.read_sbml_model("data/iJO1366.xml")
# 基本介质（同模块1）

# 示例1: 新目标 - 最大化 ATP 产量（而非生长）
atp_rxn = model.reactions.ATPM  # ATP 维护反应
original_obj = model.objective  # 保存原目标
model.objective = atp_rxn  # 设置新目标
solution_atp = model.optimize()

BIOMASS_ID = 'BIOMASS_Ec_iJO1366_core_53p95M'

print(f"原生长目标: {original_obj} -> 通量 {solution_atp.fluxes[BIOMASS_ID]:.4f}")
print(f"ATP 最大产量: {solution_atp.objective_value:.4f} mmol/gDW/h")

# 示例2: 多目标 - 生长与产物耦合（权重求和）
product_rxn = model.reactions.EX_succ_e
product_rxn.bounds = (0, 1000)
coupled_obj = 0.8 * model.objective.expression + 0.2 * product_rxn.flux_expression  # 线性组合
model.objective = coupled_obj
solution_coupled = model.optimize()

# 可视化：目标函数比较（雷达图）
categories = ['Growth', 'ATP Yield', 'Succinate']
values_original = [0.87, 0, 0]  # 简化示例值
values_atp = [0.5, 12.0, 0]  # 假设值
values_coupled = [0.7, 0, 5.0]
angles = np.linspace(0, 2*np.pi, len(categories), endpoint=False).tolist()
values_original += values_original[:1]
values_atp += values_atp[:1]
values_coupled += values_coupled[:1]
angles += angles[:1]

fig, ax = plt.subplots(figsize=(6, 6), subplot_kw=dict(projection='polar'))
ax.plot(angles, values_original, 'o-', linewidth=2, label='Original')
ax.fill(angles, values_original, alpha=0.25)
ax.plot(angles, values_atp, 'o-', linewidth=2, label='Max ATP')
ax.fill(angles, values_atp, alpha=0.25)
ax.plot(angles, values_coupled, 'o-', linewidth=2, label='Coupled')
ax.fill(angles, values_coupled, alpha=0.25)
ax.set_xticks(angles[:-1])
ax.set_xticklabels(categories)
ax.set_ylim(0, max(max(values_original), max(values_atp), max(values_coupled)) * 1.1)
plt.title('Objective Function Comparison')
plt.legend(loc='upper right', bbox_to_anchor=(0.1, 0.1))
plt.show()

# 重置目标
model.objective = original_obj




# module 9 两个变量关系分析
import cobra
import matplotlib.pyplot as plt
import numpy as np

# 加载模型
model = cobra.io.read_sbml_model("data/iJO1366.xml")  # 替换为实际路径

# 生长反应 ID
BIOMASS_ID = 'BIOMASS_Ec_iJO1366_core_53p95M'

# 固定糖耗：葡萄糖摄取 -10 mmol/gDW/h
model.reactions.EX_glc__D_e.bounds = (-10, 0)

# 必需营养（固定）
model.reactions.EX_nh4_e.bounds = (-20, 0)  # 氮源
model.reactions.EX_pi_e.bounds = (-20, 0)    # 磷源
model.reactions.EX_so4_e.bounds = (-20, 0)   # 硫源
model.reactions.EX_h_e.bounds = (-1000, 1000)  # 质子平衡
model.reactions.EX_h2o_e.bounds = (-1000, 1000)  # 水

# 调试：初始 FBA（需氧高值，预期生长 ~0.87）
model.reactions.EX_o2_e.bounds = (-30, 0)
initial_solution = model.optimize()
print(f"初始需氧生长率: {initial_solution.objective_value:.4f} h⁻¹")

# 步骤2: 变异氧气摄入，计算生长率
oxygen_levels = np.arange(0, 31, 2)  # 从 0 (厌氧) 到 30 mmol/gDW/h，步长 2（减少计算时间）
growth_rates = []

for o2_uptake in oxygen_levels:
    model.reactions.EX_o2_e.bounds = (-o2_uptake, 0)  # 设置氧气摄入（负值为摄取）
    solution = model.optimize()
    growth = solution.objective_value
    growth_rates.append(growth)
    print(f"氧气摄入 {o2_uptake:.1f} mmol/gDW/h -> 生长率: {growth:.4f} h⁻¹")

# 步骤3: 可视化
plt.figure(figsize=(8, 6))
plt.plot(oxygen_levels, growth_rates, 'o-', color='blue', linewidth=2, markersize=4)
plt.xlabel('O2 uptake rate (mmol/gDW/h)')
plt.ylabel('growth rate (h⁻¹)')
plt.grid(True, alpha=0.3)
plt.legend()
plt.tight_layout()
plt.show()

# 输出数据表（可选，便于分析）
import pandas as pd
df = pd.DataFrame({'Oxygen Uptake': oxygen_levels, 'Growth Rate': growth_rates})
print("\n数据表:")
print(df.round(4))



# new task ############

# 突变株的表型模拟 (MOMA)
# module 10 MOMA vs FBA for Mutants
import cobra
from cobra.flux_analysis import moma

print("\n=== Module 12: FBA vs. MOMA for Knockouts ===")
model = cobra.io.read_sbml_model("data/iJO1366.xml")

# 1. 计算野生型 (WT) 参照通量
wt_solution = model.optimize()

# 2. 模拟丙酮酸激酶 (Pyruvate Kinase, PYK) 敲除
# 这是一个中心碳代谢的关键反应
mutant_model = model.copy()
mutant_model.reactions.PYK.bounds = (0, 0)

# 3. 运行MOMA算法，比较基因敲除后，代谢通量差异



# module 11 pFBA vs FBA
import cobra
from cobra.flux_analysis import pfba

print("\n=== Module 11: Parsimonious FBA (pFBA) ===")
model = cobra.io.read_sbml_model("data/iJO1366.xml")

# 1. 运行普通 FBA
fba_solution = model.optimize()
fba_flux_sum = abs(fba_solution.fluxes).sum()

# 2. 运行pFBA，比较二者异同



# module 12 Gap Filling (模拟模型修复)
import cobra
from cobra.flux_analysis import find_blocked_reactions
from cobra.flux_analysis import gapfilling

print("=== Gap Filling & Model Repair ===")
model = cobra.io.read_sbml_model("data/iJO1366.xml")
original = model.copy()
# 1. 故意破坏模型：删除用于合成必需氨基酸（如精氨酸）的关键反应
# ArgG: Argininosuccinate synthase
broken_rxn_id = 'ARGSS'
print(f"正在删除关键反应: {broken_rxn_id}...")
rxn_to_remove = model.reactions.get_by_id(broken_rxn_id)
model.remove_reactions([rxn_to_remove])

# 2. 验证模型是否“死亡”
print(f"破坏后生长率: {model.optimize().objective_value:.4f}")  # 预期接近 0

# 3. 诊断：寻找阻塞反应 (Blocked Reactions)
# 注意：如果一个通路上游断了，下游很多反应都会变成 blocked
blocked = find_blocked_reactions(model)
print(f"当前被阻塞的反应数量: {len(blocked)}")

# 4. 采用gap filling, 或者添加新反应填补缺口
gap_solutions = gapfilling(
        model=model,
        universe=original,  # 候选反应集
        growth=0.1,         # 目标生长率阈值（至少 0.1 h⁻¹）
        max_results=1,      # 返回前1个最佳方案
        max_iterations=100  # 最大迭代（避免超时）
    )


