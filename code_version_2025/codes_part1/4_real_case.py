# -*- coding: utf-8 -*-
# this code is used to conduct flux sampling analysis for Streptomyces
# Hongzhong Lu
# started
from __future__ import print_function
import os
# for the sampling analysis, You have uncovered a small bug :) Until it's fixed, please downgrade to cobra 0.14.2 and remove the line referencing model._sbml
from cobra.io import read_sbml_model
from cobra.flux_analysis import pfba
from cobra.sampling.optgp import OptGPSampler
import pprint
import sys
import pandas as pd
from cobra import Model, Reaction, Metabolite


# test for the GEM_new
GEM_new = read_sbml_model('./data/iKS1317.xml')
model1 = GEM_new.copy()
model = model1

"""
GEM_new2 = read_sbml_model('../data/Sco-GEM.xml')
"""

# firstly check the cell growth
# flux calculation under reference condition
model.reactions.get_by_id("EX_glc__D_e").bounds = (-12.11,-12.11) #glu
model.reactions.get_by_id("EX_co2_e").bounds = (24.48,24.48) #co2
# optimize
model.objective = 'BIOMASS_SCO'
pfba_solution = pfba(model)
flux1 = pd.DataFrame(pfba_solution.fluxes)
flux1.index.name = 'rxnID'
flux1.reset_index(inplace=True)
# check the exchange rates
flux_exchange = flux1[flux1['rxnID'].str.contains('EX_')]
pfba_solution.fluxes['BIOMASS_SCO']


'''
# sampling
# GLPK is good there since it has faster basis recycling than Gurobi or Cplex.
gr = pfba_solution.fluxes['BIOMASS_SCO']
co = model.problem.Constraint(model.reactions.r_4041.flux_expression, lb=0.95*(gr))
model.add_cons_vars([co])
optgp = OptGPSampler(model, processes = 8)
s = optgp.sample(1000)
'''


# just for the analysis
# here we fist fix the biomass and maximize the production
model = read_sbml_model('./data/iKS1317.xml')
# Add the reactions for TAG degrade
model.add_metabolites(Metabolite('tag160_e'))
# TAG degradation
# this initial function may be not right, should check
dict1 = {model.metabolites.get_by_id('tag160_c'): -1,
        model.metabolites.get_by_id('coa_c'): -24,
        model.metabolites.get_by_id('atp_c'): -3,
        model.metabolites.get_by_id('fad_c'): -21,
        model.metabolites.get_by_id('nad_c'): -21,
        model.metabolites.get_by_id('h2o_c'): -21,
        model.metabolites.get_by_id('nadh_c'):21,
        model.metabolites.get_by_id('accoa_c'):24,
        model.metabolites.get_by_id('fadh2_c'):21,
        model.metabolites.get_by_id('glyc_c'):1,
        model.metabolites.get_by_id('amp_c'):3,
        model.metabolites.get_by_id('ppi_c'):3,
        model.metabolites.get_by_id('h_c'):24
        }

# using the general procedures to add the new function
reaction = Reaction('TAG_degrade')
reaction.name = 'TAG_degradation'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict1)
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

# transport
dict2 = {model.metabolites.get_by_id('tag160_e'): -1,
        model.metabolites.get_by_id('tag160_c'):1,
        }
reaction = Reaction('TAG_transport')
reaction.name = 'TAG_transport'
reaction.subsystem = 'new added'
reaction.lower_bound = 0.  # This is the default
reaction.upper_bound = 1000.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict2)
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])

# exchange
dict3 = {model.metabolites.get_by_id('tag160_e'): -1
        }
reaction = Reaction('TAG_exchange')
reaction.name = 'TAG_exchange'
reaction.subsystem = 'new added'
reaction.lower_bound = -1000.  # This is the default
reaction.upper_bound = 0.  # This is the default
reaction.EC = ''
reaction.add_metabolites(dict3)
reaction.gene_reaction_rule = ''
model.add_reactions([reaction])
# check the added reactions
for r in model.reactions:
    print(r.id, r.name, r.reaction, r.compartments, sep="\t")


# with glucose uptake as objective function
# calculation
#M145 no feeding
with model:
    model.reactions.get_by_id("EX_glc__D_e").bounds = (-1000, 0)  # glu
    model.reactions.get_by_id("TAG_exchange").bounds = (-0.676, -0.676)  # fix TAG
    model.reactions.get_by_id('BIOMASS_SCO').bounds = (1, 1)
    model.reactions.get_by_id("EX_co2_e").bounds = (24.48, 24.48)  # co2
    model.reactions.get_by_id("DM_ACT_c").bounds = (1.1425, 1.1425)

    #model.objective = 'DM_ACT_c'
    model.objective = 'EX_glc__D_e'
    pfba_solution = pfba(model)
    #print('act:', pfba_solution.fluxes['DM_ACT_c'])
    print('glc:', pfba_solution.fluxes['EX_glc__D_e'])
    #flux sampling
    gr1 = pfba_solution.fluxes['EX_glc__D_e']
    co = model.problem.Constraint(model.reactions.EX_glc__D_e.flux_expression, lb=1.02 * (gr1))
    model.add_cons_vars([co])
    optgp = OptGPSampler(model, processes=8)
    s1 = optgp.sample(1000)



#HY01 no feeding
with model:
    model.reactions.get_by_id("EX_glc__D_e").bounds = (-1000, 0)  # glu
    model.reactions.get_by_id("TAG_exchange").bounds = (-1.438, -1.438)  # fix TAG
    model.reactions.get_by_id('BIOMASS_SCO').bounds = (1, 1)
    model.reactions.get_by_id("EX_co2_e").bounds = (13.51, 13.51)  # co2
    model.reactions.get_by_id("DM_ACT_c").bounds = (2.579, 2.579)
    #model.objective = 'DM_ACT_c'
    model.objective = 'EX_glc__D_e'
    pfba_solution = pfba(model)
    #print('act:', pfba_solution.fluxes['DM_ACT_c'])
    print('glc:', pfba_solution.fluxes['EX_glc__D_e'])
    #flux sampling
    gr2 = pfba_solution.fluxes['EX_glc__D_e']
    co = model.problem.Constraint(model.reactions.EX_glc__D_e.flux_expression, lb=1.02 * (gr2))
    model.add_cons_vars([co])
    optgp = OptGPSampler(model, processes=8)
    s2 = optgp.sample(1000)


import seaborn as sns
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(4,3))
sns.distplot(s1['CS'])
sns.distplot(s2['CS'])
fig.legend(labels=['M145, no feeding','HY01, no feeding'])
plt.show()

fig = plt.figure(figsize=(4,3)) # when re-run this step, the pycharm is stopped unexpectedly
sns.distplot(s1['ACONTa'])
sns.distplot(s2['ACONTa'])
fig.legend(labels=['M145, no feeding','HY01, no feeding'])
plt.show()



#M145 feeding
with model:
    model.reactions.get_by_id("EX_glc__D_e").bounds = (-1000, 0)  # glu
    model.reactions.get_by_id("TAG_exchange").bounds = (-0.476, -0.476)  # fix TAG
    model.reactions.get_by_id('BIOMASS_SCO').bounds = (8.4, 8.4)
    model.reactions.get_by_id("EX_co2_e").bounds = (261.16, 261.16)  # co2
    model.reactions.get_by_id("DM_ACT_c").bounds = (1.459, 1.459)

    #model.objective = 'DM_ACT_c'
    model.objective = 'EX_glc__D_e'
    pfba_solution = pfba(model)
    #print('act:', pfba_solution.fluxes['DM_ACT_c'])
    print('glc:', pfba_solution.fluxes['EX_glc__D_e'])
    #flux sampling
    gr1 = pfba_solution.fluxes['EX_glc__D_e']
    co = model.problem.Constraint(model.reactions.EX_glc__D_e.flux_expression, lb=1.02 * (gr1))
    model.add_cons_vars([co])
    optgp = OptGPSampler(model, processes=8)
    s3 = optgp.sample(1000)



#HY01 feeding
with model:
    model.reactions.get_by_id("EX_glc__D_e").bounds = (-1000, 0)  # glu
    model.reactions.get_by_id("TAG_exchange").bounds = (-1.088, -1.088)  # fix TAG
    model.reactions.get_by_id('BIOMASS_SCO').bounds = (8.3, 8.3)
    model.reactions.get_by_id("EX_co2_e").bounds = (205.31, 205.31)  # co2
    model.reactions.get_by_id("DM_ACT_c").bounds = (3.5125, 3.5125)
    #model.objective = 'DM_ACT_c'
    model.objective = 'EX_glc__D_e'
    pfba_solution = pfba(model)
    #print('act:', pfba_solution.fluxes['DM_ACT_c'])
    print('glc:', pfba_solution.fluxes['EX_glc__D_e'])
    #flux sampling
    gr2 = pfba_solution.fluxes['EX_glc__D_e']
    co = model.problem.Constraint(model.reactions.EX_glc__D_e.flux_expression, lb=1.02 * (gr2))
    model.add_cons_vars([co])
    optgp = OptGPSampler(model, processes=8)
    s4 = optgp.sample(1000)


import seaborn as sns
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(4,3))
sns.distplot(s3['CS'])
sns.distplot(s4['CS'])
fig.legend(labels=['M145, feeding','HY01, feeding'])
plt.show()

fig = plt.figure(figsize=(4,3))
sns.distplot(s3['ACONTa'])
sns.distplot(s4['ACONTa'])
fig.legend(labels=['M145, feeding','HY01, feeding'])
plt.show()

