# -*- coding: utf-8 -*-
'''Use ecFSEOF algorithm with ecGEM to predict potential targets for spermidine production.'''
import cobra
import sys

# add the path of the codebase
sys.path.append('codebase')
import codebase


# 1.load model
model=cobra.io.read_sbml_model('model/ecYeastGEM_batch.xml')

# 2. Set the parameters
model_type='ecGEM'
growth='r_2111'    # biomass reaction
product='r_2051'   # spermidine exchange reaction
c_source = 'r_1714_REV'   # glucose exchange reaction
c_uptake = 1

# 3. Set the scanning range based on the theoretical maximum product yield
with model:
    model.objective=growth  # biomass rxn
    model.reactions.get_by_id(c_source).bounds=c_uptake,c_uptake
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
final_result.to_csv('result/ecfseof_output.csv')
