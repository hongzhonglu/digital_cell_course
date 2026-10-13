# -*- coding: utf-8 -*-
'''Use optForce algorithm with GEM to predict potential targets for spermidine production.'''
import cobra
import sys

# add the path of the codebase
sys.path.append('codebase')
import codebase


# 1.load model
model=cobra.io.read_sbml_model('model/yeast-GEM.xml')


# 2. Set the parameters
model_type='GEM'
growth='r_2111'
product='r_2051'   # spermidine exchange reaction
c_source = 'r_1714'
c_uptake = 1 #mmol/gDW.h

biomass_const = 0.9
target_const = 0.9

bio_max=model.slim_optimize() #  /h


# 3. Simulation
# for wild type
print('FVA for wild type')
bio_fva_result = codebase.run_FVA(model, target_rxn=growth, target_const_rate=biomass_const)

# for production condition
flux_constraint = {growth:(0.1*bio_max, 0.1*bio_max)}
print('FVA for production condition')
target_fva_result = codebase.run_FVA(model, target_rxn=product, target_const_rate=target_const, flux_constraints=flux_constraint)

down_targets, up_targets = codebase.compare_FV_range(bio_fva_result, target_fva_result)

# 4. Get gene targets result
final_result=codebase.parse_optforce_result(model,down_rxns=down_targets,
                                            up_rxns=up_targets)

print(f'Total predicted targets with {model_type} by optForce:{len(final_result)}')
print(f'Overexpression targets:{len(final_result.loc[final_result["action"] == "OE"])}')
print(f'Knockdown targets:{len(final_result.loc[final_result["action"] == "KD"])}')

# save result
final_result.to_csv(r'result/optForce_output.csv')



