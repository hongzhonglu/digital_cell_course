# 导入包
import cobra  #导入cobra包
import pandas as pd



# 导入模型1
file = 'data/yeast-GEM.xml'  #你所需要导入模型地址
file2 = 'data/yeast-GEM.mat'  #你所需要导入模型地址
model = cobra.io.read_sbml_model(file)  #读取模型
model = cobra.io.load_matlab_model(file2)



# 导入模型2
file = 'data/e_coli_core.xml'  #你所需要导入模型地址
ecoli = cobra.io.read_sbml_model(file)  #读取模型

# 导入模型3
file = 'data/mini_cobra.xml'  #你所需要导入模型地址
mini = cobra.io.read_sbml_model(file)  #读取模型




# 查看模型
print(len(model.reactions))#输出模型中反应数量
print(len(model.metabolites))#输出模型中代谢物数量
print(len(model.genes))#输出模型中基因数量
print(model.reactions.r_0001.annotation)#反应注释
print(model.reactions.r_0001.name)#反应name
print(model.reactions.r_0001.id)#反应id
print(model.reactions.r_0001.gene_reaction_rule)#反应基因注释
print(model.reactions.r_0001.reaction)#反应内容
print(model.reactions.r_0001.upper_bound)#反应上界
print(model.reactions.r_0001.lower_bound)#反应下界



#简单分析
model.objective='r_2111'#设置模型目标反应
fba=model.optimize()#进行FBA模拟分析
print(fba.objective_value)#输出求解的目标反应的通量
print(fba.fluxes)#输出求解的所有的反应的通量
model.reactions.query("biomass",'name')#找到特定name反应
fba_frame = fba.to_frame()
print(fba_frame[fba_frame.fluxes.abs() > 1e-4])#输出所有通量大于0.0001 mmol/g*DWh-1的反应



#模拟敲除
Target_gene='Q0045'#模拟敲除的基因
with model:
    model.genes.get_by_id(Target_gene).knock_out()#模拟敲除基因
    fba_simulating_with_deletion = model.optimize()#模拟敲除后进行FBA模拟分析
    print(fba_simulating_with_deletion.objective_value)#输出求解目标反应的通量



# 约束模型并求解-FBA & pFBA
model.reactions.get_by_id("r_1714").bounds = (-3.255,-3.255) #glu
model.reactions.get_by_id("r_1808").bounds = (0, 0) #glyc
model.reactions.get_by_id("r_2033").bounds = (0.047, 0.047) #pyr
model.reactions.get_by_id("r_2056").bounds = (0,0) #succ
model.reactions.get_by_id("r_1634").bounds = (1.397,1.397) #ac
model.reactions.get_by_id("r_1697").bounds = (7.633,7.633) #co2
# optimize
solution = model.optimize()
model.summary()
fluxes = pd.DataFrame(solution.fluxes)

pfba_solution = cobra.flux_analysis.pfba(model, fraction_of_optimum=0.95)#pFBA模拟分析，fraction_of_optimum参数为pFBA的objective反应通量与FBA的objective反应通量的的比例
print(solution.fluxes.r_2111)
print(solution.fluxes.abs().sum())#FBA分析的反应的通量绝对值总量
print(pfba_solution.fluxes.r_2111)
print(pfba_solution.fluxes.abs().sum())#pFBA分析的反应通量绝对值总量（pFBA最小化该值）
