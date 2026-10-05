############################################
#  ecYeast simulation constrained with total protein pool
#  Hongzhong Lu
#  2019-7-15
#############################################
from cobra.io import load_matlab_model, read_sbml_model
import os
import sys


# import self function
from in_silico_strain_design_scripts import *


# we further try to use the general yeastGEM model to do the above analysis
# load the GEM model
ecYeast = load_matlab_model("../data/wMC1.mat")
ecYeast.solver ='glpk' # if not setting as this, the model correction will produce error
model1 = ecYeast.copy()

# initial optimize
model = ecYeastMinimalMedia(model1)

# time point at 20.52 for the wild
# fix the growth rate
add_constraint = dict()
add_constraint['r_1761_REV'] = 0.07763 # ethanol secretion
add_constraint['r_2111'] = 0.049252991 # growth rate
for key in add_constraint:
    print(key)
    model.reactions.get_by_id(key).bounds = (add_constraint[key], add_constraint[key])

model.reactions.get_by_id('r_1714_REV').bounds = (0, 1000) # open the glucose uptake rate

# set the glucose uptake as the objective funcition: minimaztion
model.objective = {model.reactions.get_by_id('r_1714_REV'): -1}
# second optimize
solution2 = model.optimize()
model.summary()
# then fix the glucose uptake rate
qs = solution2.fluxes['r_1714_REV']
model.reactions.get_by_id("r_1714_REV").bounds = (0.99*qs, qs) # open the ethanol uptake rate
# then set the protein pool as the objective function, minimizing it
model.objective = {model.reactions.prot_pool_exchange: -1}
solution3 = model.optimize()
# check the result
solution3.fluxes['r_2111'] # biomass
solution3.fluxes['r_1992_REV'] # O2
solution3.fluxes['r_1672'] # co2
solution3.fluxes['prot_pool_exchange'] # protein pool



