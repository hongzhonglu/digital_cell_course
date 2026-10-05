###################################
#  Simulation with yeast GEM
#  Hongzhong Lu
#  2019-7-15
###################################
from cobra.io import load_matlab_model, read_sbml_model
import os
import sys

sys.path.append('/Users/xluhon/Google Drive/AA. 研究生课件准备/AA.课件汇总/v2 part2/10_code_practice/code/') # 找到绝对路径
os.chdir('/Users/xluhon/Google Drive/AA. 研究生课件准备/AA.课件汇总/v2 part2/10_code_practice/code/') # 找到绝对路径





# import self function
from mainFunction import *
from in_silico_strain_design_scripts import *


# we further try to use the general yeastGEM model to do the above analysis
# load the GEM model
yeastGEM = read_sbml_model("../data/yeastGEM.xml")
yeastGEM.solver ='glpk' # if not setting as this, the model correction will produce error
model = yeastGEM.copy()
# Correct metabolite ids:
model.id = model.id.replace('__46__', '_')
for met in model.metabolites:
    met.id = met.id.replace('__91__','_')
    met.id = met.id.replace('__93__','')



# Add the constraint
# set the bounds
# time point at 20.52 for the wild
model.reactions.get_by_id("r_1714").bounds = (-0.49,-0.49) #glu
model.reactions.get_by_id("r_1761").bounds = (-0.07763381,-0.07763381) #ethonal
# optimize
solution = model.optimize()
model.summary()




# obtain the oxygen and CO2
solution.fluxes['r_2111'] # biomass
solution.fluxes['r_1992'] # O2
solution.fluxes['r_1672'] # co2



# time point at 42.4 for the wild
model.reactions.get_by_id("r_1714").bounds = (-0.503,-0.503) #glu
model.reactions.get_by_id("r_1761").bounds = (-0.0794754,-0.0794754) #ethonal

# optimize
solution = model.optimize()
model.summary()
# obtain the oxygen and CO2
solution.fluxes['r_2111'] # biomass
solution.fluxes['r_1992'] # O2
solution.fluxes['r_1672'] # co2




'''
a small test
# xylose as the carbon source
model.reactions.get_by_id("r_1714").bounds = (0,0) #glu
model.reactions.get_by_id("r_1718").bounds = (-5, -5) #ethonal
# optimize
solution = model.optimize()
model.summary()

model_no_ADH6 = getModelWithRemoveGene(model0= model, gene_remove0=['YMR318C'])
solution = model.optimize()
model.summary()

model_no_GPD1 = getModelWithRemoveGene(model0= model, gene_remove0=['YDL022W'])
solution = model.optimize()
model.summary()

model_no_ALD4 = getModelWithRemoveGene(model0= model, gene_remove0=['YOR374W'])
solution = model.optimize()
model.summary()
'''
