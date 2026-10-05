# -*- coding: utf-8 -*-
'''Use FSEOF algorithm with GEM to predict potential targets for spermidine production.'''
import cobra
import sys

# add the path of the codebase
sys.path.append('codebase')
import codebase


# 1.load model
model=cobra.io.read_sbml_model('model/yeast-GEM.xml')

itaconate_c = cobra.Metabolite('s_itaconate_c', name='Itaconate[c]', compartment='c')
itaconate_e = cobra.Metabolite('s_itaconate_e', name='Itaconate[e]', compartment='e')

cis_aconitate = model.metabolites.get_by_id("s_0529")
proton = model.metabolites.get_by_id("s_0794")
co2 = model.metabolites.get_by_id("s_0450")

# cis-aconitate + H+ -> CO2 + Itaconate
r_cad = cobra.Reaction('R_CAD')
r_cad.name = 'Cis-aconitic acid decarboxylase'
r_cad.add_metabolites({
    cis_aconitate: -1,
    proton: -1,
    co2: 1,
    itaconate_c: 1
})

# --- 反应 2: 转运 (胞内 -> 胞外) ---
r_trans = cobra.Reaction('R_ITA_tr')
r_trans.name = 'Itaconic acid transport'
r_trans.add_metabolites({
    itaconate_c: -1,
    itaconate_e: 1
})

# --- 反应 3: 交换 (排出到环境) ---
r_ex = cobra.Reaction('EX_itaconate')
r_ex.name = 'Itaconic acid exchange'
r_ex.add_metabolites({
    itaconate_e: -1
})

# 将这三个反应加入模型
model.add_reactions([r_cad, r_trans, r_ex])

# 5. 设定目标并计算
model.objective = 'EX_itaconate'  # 将目标函数设为衣康酸出口
solution = model.optimize()

# 6. 输出结果
print("-" * 30)
print(f"衣康酸最大合成速率: {solution.objective_value:.4f} mmol/gDW/h")
print("-" * 30)





# 2. Set the parameters
model_type='GEM'
growth='r_2111'
product= 'EX_itaconate' #'r_2051'   # spermidine exchange reaction
c_source = 'r_1714'
c_uptake = 1  # mmol/gDW.h

# 3. Set the scanning range based on the theoretical maximum product yield
with model:
    model.objective=growth  # biomass rxn
    model.reactions.get_by_id(c_source).bounds=-c_uptake,-c_uptake
    gluc_MW=0.180156  # g/mmol
    max_yield=model.slim_optimize()/(c_uptake*gluc_MW) # gDW / gGluc
    expYield=max_yield*0.49
    alphaLims=(0.5*expYield,2*expYield)  # The scanning range of the product yield can be adjusted according to the actual situation


# set the action thresholds
action_thresholds=[0.05,0.3,1.05]     # rules for knockout, knockdown and overexpression

# 4. Run the FSEOF algorithm to identify the potential targets
Nsteps = 16  # number of FBA steps in ecFSEOF
results = codebase.run_FSEOF(model=model,
                          targetID=product,
                          c_source=c_source,
                          c_uptake=c_uptake,
                          alphaLims=alphaLims,
                          Nsteps=Nsteps,
                          model_type=model_type)

# Format results table
final_result = results['geneTable']
final_result.loc[final_result['k_score'] >= action_thresholds[2], 'action'] = 'OE'
final_result.loc[final_result['k_score'] <= action_thresholds[1], 'action'] = 'KD'
final_result.loc[final_result['k_score'] <= action_thresholds[0], 'action'] = 'KO'
# remove genes with no action
# final_result = final_result.loc[final_result['action'].notnull()]
final_result = final_result.loc[final_result['action'].isin(['OE', 'KD', 'KO'])]
print(f'Total predicted targets with {model_type}:{len(final_result)}')
print(f'Overexpression targets:{len(final_result.loc[final_result["action"] == "OE"])}')
print(f'Knockdown targets:{len(final_result.loc[final_result["action"] == "KD"])}')
print(f'Knockout targets:{len(final_result.loc[final_result["action"] == "KO"])}')

# save result
final_result.to_csv('result/fseof_output.csv')
