###################################################
#  Add carotenoid synthesis reaction into ecYeast
#  Hongzhong Lu
#  2019-7-15
###################################################
from cobra.io import load_matlab_model, read_sbml_model, write_sbml_model # when using save_matlab_model, the compartment information in the new model will be lost!!!
from cobra import Reaction, Metabolite

import os
import sys
import pandas as pd


# import self function
from mainFunction import *
from in_silico_strain_design_scripts import *


# we further try to use the general yeastGEM model to do the above analysis
# load the GEM model

# an example
model_name = 'wProt1'

model_input = '../data/' + model_name + '.mat'
model_output = '../data/' + model_name + '_new'+'.xml'
ecModel = load_matlab_model(model_input)
model = ecModel.copy()
# check the metabolite name
for met in model.metabolites:
    print(met.id)

# kcat calculation for  crtYB
# EC 2.5.1.32
MW1 = 74736
specific_activity1 =0.067
kcat1 = specific_activity1/60*MW1/1000

# EC 5.5.1.19
MW2 = 74736
specific_activity2 = 0.356  # unit:µmol/min/mg
kcat2 = specific_activity2/60*MW2/1000 # unit: /s




# firstly add the new metabolite
model.add_metabolites(Metabolite('s_4207', compartment='c', formula= 'C40H64', name='15-cis-phytoene'))
model.add_metabolites(Metabolite('s_4208',compartment='c', formula= 'C40H56', name='all-trans-lycopene'))
model.add_metabolites(Metabolite('s_4209', compartment='c', formula= 'C40H56', name='γ-carotene'))
model.add_metabolites(Metabolite('s_4210',compartment='c', formula= 'C40H56', name='β-carotene'))
# add a new protein into the model
model.add_metabolites(Metabolite('prot_crtYB',compartment='c', formula= '', name=''))
model.add_metabolites(Metabolite('product_all',compartment='c', formula= '', name=''))




# define the reactions
dict1 = {model.metabolites.get_by_id('s_0189'): -2,
         model.metabolites.get_by_id('prot_crtYB'): -1/(kcat1*3600),
         model.metabolites.get_by_id('s_4207'):1,
         model.metabolites.get_by_id('s_0633'):2
        }


dict2 = {model.metabolites.get_by_id('s_4207'): -1,
         model.metabolites.get_by_id('s_1207'):-4,
         model.metabolites.get_by_id('s_4208'):1,
         model.metabolites.get_by_id('s_1212'):4
        }


dict3 = {model.metabolites.get_by_id('s_4208'): -1,
         model.metabolites.get_by_id('prot_crtYB'): -1/(kcat2*3600),
         model.metabolites.get_by_id('s_4209'):1
        }

dict4 = {model.metabolites.get_by_id('s_4209'): -1,
         model.metabolites.get_by_id('prot_crtYB'): -1/(kcat2*3600),
         model.metabolites.get_by_id('s_4210'):1
        }


# exchange
dict5 = {model.metabolites.get_by_id('s_4207'): -1
        }

dict6 = {model.metabolites.get_by_id('s_4208'): -1
        }

dict7 = {model.metabolites.get_by_id('s_4210'): -1
        }

# add the protein pool to the crtYB
dict8 = {model.metabolites.get_by_id('prot_pool'): -MW1/1000,
         model.metabolites.get_by_id('prot_crtYB'):1
        }


# add the sum of product
dict9 = {model.metabolites.get_by_id('s_4207'): -1,
        model.metabolites.get_by_id('s_4208'): -1,
        model.metabolites.get_by_id('s_4210'): -1,
        model.metabolites.get_by_id('product_all'): 1
        }
dict10 = {model.metabolites.get_by_id('product_all'): -1
        }



# using the general procedures to add the new function
# using a loop to add all the reactions into yeast-GEM
# input the original data
descriptionRxn = pd.read_excel('../data/dd-decaf dsm carotenoid data heterologous reactions.xlsx', 'rxn_ec')
compositionRxn = [dict1, dict2, dict3, dict4, dict5, dict6, dict7,dict8, dict9,dict10]

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
    #if row['GPR'] is not None:
    #    reaction.gene_reaction_rule = row['GPR']
    #else:
    #    reaction.gene_reaction_rule = ''

    if row['GPR'] =="None":
        reaction.gene_reaction_rule = row['GPR']
    else:
        reaction.gene_reaction_rule = ''
    model.add_reactions([reaction])

# check the reaction
for r in model.reactions:
    print(r.id, r.name, r.reaction, r.compartments, sep="\t")

# save the model
write_sbml_model(model, model_output)




