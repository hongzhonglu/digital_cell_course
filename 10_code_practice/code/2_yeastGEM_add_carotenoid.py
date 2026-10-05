###################################################
#  Add carotenoid synthesis reaction into yeast GEM
#  Hongzhong Lu
#  2019-7-15
###################################################
from cobra.io import load_matlab_model, read_sbml_model, write_sbml_model
from cobra import  Reaction, Metabolite

import os
import sys
import pandas as pd


# import self function
from mainFunction import *
from in_silico_strain_design_scripts import *


# we further try to use the general yeastGEM model to do the above analysis
# load the GEM model
yeastGEM = read_sbml_model("../data/yeastGEM.xml")
yeastGEM.solver ='glpk' # if not setting as this, the model correction will produce error
model = yeastGEM.copy()
# Correct metabolite ids:
# note: this replace is not used any more
model.id = model.id.replace('__46__', '_')
for met in model.metabolites:
    met.id = met.id.replace('__91__','_')
    met.id = met.id.replace('__93__','')
    print(met.id)

# firstly add the new metabolite
model.add_metabolites(Metabolite('s_4207[c]', compartment='c', formula= 'C40H64', name='15-cis-phytoene'))
model.add_metabolites(Metabolite('s_4208[c]',compartment='c', formula= 'C40H56', name='all-trans-lycopene'))
model.add_metabolites(Metabolite('s_4209[c]', compartment='c', formula= 'C40H56', name='γ-carotene'))
model.add_metabolites(Metabolite('s_4210[c]',compartment='c', formula= 'C40H56', name='β-carotene'))


# define the reactions
# reaction 1
dict1 = {model.metabolites.get_by_id('s_0189[c]'): -2,
         model.metabolites.get_by_id('s_4207[c]'):1,
         model.metabolites.get_by_id('s_0633[c]'):2
        }


dict2 = {model.metabolites.get_by_id('s_4207[c]'): -1,
         model.metabolites.get_by_id('s_1207[c]'):-4,
         model.metabolites.get_by_id('s_4208[c]'):1,
         model.metabolites.get_by_id('s_1212[c]'):4
        }


dict3 = {model.metabolites.get_by_id('s_4208[c]'): -1,
         model.metabolites.get_by_id('s_4209[c]'):1
        }

dict4 = {model.metabolites.get_by_id('s_4209[c]'): -1,
         model.metabolites.get_by_id('s_4210[c]'):1
        }


# exchange
dict5 = {model.metabolites.get_by_id('s_4207[c]'): -1
        }

dict6 = {model.metabolites.get_by_id('s_4208[c]'): -1
        }

dict7 = {model.metabolites.get_by_id('s_4210[c]'): -1
        }


# using the general procedures to add the new function
# using a loop to add all the reactions into yeast-GEM
# input the original data
descriptionRxn = pd.read_excel('../data/dd-decaf dsm carotenoid data heterologous reactions.xlsx', 'rxn')
descriptionRxn = descriptionRxn[~descriptionRxn["rxn"].isna()]
compositionRxn = [dict1, dict2, dict3, dict4, dict5, dict6, dict7]




for index, row in descriptionRxn.iterrows():
    print(index)
    print(row['rxn'])
    reaction = Reaction(row['rxn'])
    reaction.name = ''
    reaction.subsystem = row['subsystem']
    reaction.lower_bound = row['lb']  # This is the default
    reaction.upper_bound = row['ub']  # This is the default
    reaction.EC = row['EC']
    reaction.add_metabolites(compositionRxn[index])
    if row['GPR'] =="None":
        reaction.gene_reaction_rule = row['GPR']
    else:
        reaction.gene_reaction_rule = ''
    model.add_reactions([reaction])


# check the reaction
for rxn in model.reactions:
    print(rxn.id)
    print(rxn.reaction)

# save the model
write_sbml_model(model, "../data/yeastGEM_add_carotenoid.xml")

