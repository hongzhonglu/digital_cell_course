###################################
#  application for ecYeast
#  explore the relation between qs, growth and protein pool using ecModel
###################################
from cobra.io import load_matlab_model, read_sbml_model, write_sbml_model
import os
import sys
import matplotlib.pyplot as plt # For python 3.6: sudo apt-get install python3.6-tk   TO INSTALL Tkinter


# import self function
from in_silico_strain_design_scripts import *


# load the ecYeast model
# ecYeast = load_matlab_model("../data/ec3HP.mat")
# write_sbml_model(ecYeast, "../data/ec3HP.xml")
ecYeast = read_sbml_model("../data/ec3HP.xml")


"""general analysis"""
# Open the glucose uptake
glucose_uptake = list(frange(0, 15, 0.2))
growth = []
protein_pool = []
model0 = ecYeast.copy()
model0 = ecYeastMinimalMedia(model0)
model0.objective = 'r_2111'
model0.reactions.get_by_id("r_2111").lower_bound = 0

for i in glucose_uptake:
    model0.reactions.get_by_id("r_1714_REV").bounds = (0, i) #D-glucose exchange (reversible)
    solution = model0.optimize()
    print(solution.objective_value)
    growth.append(solution.objective_value)
    protein_pool.append(solution.fluxes['prot_pool_exchange'])
# plot the relation between the glucose uptake rate and growth rate
plt.figure()
plt.plot(glucose_uptake, growth, 'ko')
plt.xlabel('Glucose uptake(mmol/gDW.h)')
plt.ylabel('Growth rate (/h)')
plt.show()

# plot the relation between the growth rate and pro_pool_exchange
plt.figure()
plt.plot(protein_pool, growth, 'ko')
plt.xlabel('protein_pool (g/gDW)') # update
plt.ylabel('Growth rate (/h)')
plt.show()



# if we fix the glucose uptake rate as 10 mmol/gDW.h
model1 = ecYeast.copy()
model1.reactions.get_by_id("r_1714_REV").bounds = (10, 10)  # D-glucose exchange (reversible)
growth = list(frange(0, 0.37, 0.01))
production = []
for i in growth:
    model1.reactions.get_by_id("r_2111").bounds = (i, i) #D-glucose exchange (reversible)
    solution = model1.optimize()
    print(solution.objective_value)
    production.append(solution.objective_value)
# plot the result
plt.figure()
plt.plot(growth, production, 'ko')
plt.xlabel('growth(/h)')
plt.ylabel('production(mmol/gDW.h)')
plt.show()



