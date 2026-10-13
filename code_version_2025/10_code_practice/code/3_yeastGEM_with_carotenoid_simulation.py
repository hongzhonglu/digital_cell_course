###################################################
#  Simulation with yeast GEM
#  Hongzhong Lu
#  2019-7-15
###################################################
from cobra.io import load_matlab_model, read_sbml_model, write_sbml_model
import os
import sys


# import self function
from mainFunction import *
from in_silico_strain_design_scripts import *


# we further try to use the general yeastGEM model to do the above analysis
# load the GEM model
yeastGEM = read_sbml_model("../data/yeastGEM_add_carotenoid.xml")
yeastGEM.solver ='glpk' # if not setting as this, the model correction will produce error
model = yeastGEM.copy()
for x in model.reactions:
    print(x.id)

# Add the constraint
# time point at 20.61 for car002_exp_feed
model.reactions.get_by_id("r_1714").bounds = (-0.521001769,-0.521001769) #glu
model.reactions.get_by_id("r_1761").bounds = (-0.082261973,-0.082261973) #ethonal
model.reactions.get_by_id("DM_phytoene").bounds = (0.000577032, 0.000577032)
model.reactions.get_by_id("DM_lycopene").bounds = (0.006405934, 0.006405934)
model.reactions.get_by_id("DM_b_carotene").bounds = (0.001123664, 0.001123664)

# optimize
solution = model.optimize()
model.summary()

# obtain the oxygen and CO2
solution.fluxes['r_2111'] # biomass
solution.fluxes['r_1992'] # O2
solution.fluxes['r_1672'] # co2


# time point at 20.61 for wild strain
model.reactions.get_by_id("r_1714").bounds = (-0.491680795, -0.491680795) #glu
model.reactions.get_by_id("r_1761").bounds = (-0.07763381, -0.07763381) #ethonal
model.reactions.get_by_id("DM_phytoene").bounds = (0, 0)
model.reactions.get_by_id("DM_lycopene").bounds = (0, 0)
model.reactions.get_by_id("DM_b_carotene").bounds = (0, 0)

# optimize
solution = model.optimize()
model.summary()

# obtain the oxygen and CO2
solution.fluxes['r_2111'] # biomass
solution.fluxes['r_1992'] # O2
solution.fluxes['r_1672'] # co2





# time point at 20.61 for car020_exp_feed
model.reactions.get_by_id("r_1714").bounds = (-0.54075,-0.54075) #glu
model.reactions.get_by_id("r_1761").bounds = (-0.08538,-0.08538) #ethonal
model.reactions.get_by_id("DM_phytoene").bounds = (0.000866625,0.000866625)
model.reactions.get_by_id("DM_lycopene").bounds = (0.006234562, 0.006234562)
model.reactions.get_by_id("DM_b_carotene").bounds = (0.001314757, 0.001314757)

# optimize
solution = model.optimize()
model.summary()

# obtain the oxygen and CO2
solution.fluxes['r_2111'] # biomass
solution.fluxes['r_1992'] # O2
solution.fluxes['r_1672'] # co2



# time point at 20.67 for car021_exp_feed
model.reactions.get_by_id("r_1714").bounds = (-0.457,-0.457) #glu
model.reactions.get_by_id("r_1761").bounds = (-0.0722,-0.0722) #ethonal
model.reactions.get_by_id("DM_phytoene").bounds = (0.0026,0.0026)
model.reactions.get_by_id("DM_lycopene").bounds = (0.00316,0.00316)
model.reactions.get_by_id("DM_b_carotene").bounds = (0.00273,0.00273)

# optimize
solution = model.optimize()
model.summary()

# obtain the oxygen and CO2
solution.fluxes['r_2111'] # biomass
solution.fluxes['r_1992'] # O2
solution.fluxes['r_1672'] # co2
