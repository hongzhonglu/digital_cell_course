# this module is mainly for ecModel simulation
# for the ecModels with organelle constraint, can it still predict Crabtree effect? In principle it should.


import sys
import os
from cobra.io import load_matlab_model, read_sbml_model
import matplotlib.pyplot as plt
import seaborn as sns
import numpy as np
import pandas as pd

# import self function
def getRxnByReactionName(model, name):
    """
    This function is used to extract the rxn id based on rxn name
    It is suitable for ecGEMs as multiple reations could use the same name

    :param model:
    :param name:
    :return:
    """
    s = []
    for rxn in model.reactions:
        if name == rxn.name:
            #print(rxn.id)
            s.append(rxn.id)
    return s
def singleMapping (description, item1, item2, dataframe=True):
    """get the single description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    # used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    index = [None]*len(item2)
    result = [None]*len(item2)
    tt = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index[i] = item1.index(item2[i])
            result[i] = description[index[i]]
        else:
            index[i] = None
            result[i] = None
    return result
def DLecModelSimulate(model, dilution_rate):
   """
   This function is used to do simulation with ecModels using kcat value from deep learning.
   :param model: a ecModel
   :param dilution_rate: a dilution rate in range of 0-0.42 /h

   :return: solution_f: fluxes datasets

   """

   dilutionrate = dilution_rate
   ecYeast = model
   if dilutionrate >= 0.4:
       ecYeast.reactions.get_by_id("EX_protein_pool").bounds = (-167.27 * dilutionrate / 0.4, 0)  # this value is further rescaled by maximal growth rate at 0.42.
   else:
       ecYeast.reactions.get_by_id("EX_protein_pool").bounds = (-167.27, 0)  # -230/0.55*0.4, this is rescaled by maximal growth rate.

   # refer to bioRxiv
   ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
              'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange',
              'oxygen exchange', 'EX_protein_pool']
   # find the related rxnID
   idx = []
   for name0 in ex_mets:
       #print(name0)
       s = getRxnByReactionName(model=ecYeast, name=name0)
       if len(s) > 1:
           print("need check")
       elif len(s) == 1:
           idx.append(s[0])

   model_tmp = ecYeast #.copy() # this is adjusted for the newly added constraints from organelle
   model_tmp.reactions.get_by_id(idx[1]).lower_bound = -1000  # glucose uptake
   model_tmp.reactions.get_by_id(idx[0]).lower_bound = dilutionrate
   model_tmp.objective = {model_tmp.reactions.r_1714: 1}  # minimize the uptake of glucose
   solution2 = model_tmp.optimize()
   # then fix glucose uptake and minimize the protein pool
   model_tmp.reactions.get_by_id(idx[1]).lower_bound = solution2.objective_value * 1.00001
   print('Glucose uptake rate: ',solution2.objective_value)
   model_tmp.reactions.get_by_id(idx[9]).lower_bound = -1000 # protein pool
   model_tmp.objective = {model_tmp.reactions.EX_protein_pool: 1}  # minimize the usage of protein pools
   solution_f = model_tmp.optimize()
   solution_f.fluxes["EX_protein_pool"]
   solution_f.fluxes["r_1714"]
   return solution_f



# using the manual curated ecYeast from deep learning
# in this version of model, we curate the kcat for some enzymes
ecYeast = read_sbml_model("../data/ecYeast_DL_update_some_kcat.xml")






# simulate crabtree effect!
# refer to bioRxiv
ex_mets = ['biomass pseudoreaction', 'D-glucose exchange', 'acetate exchange', 'ethanol exchange',
           'glycerol exchange', 'pyruvate exchange', 'ethyl acetate exchange', 'carbon dioxide exchange', 'oxygen exchange', 'EX_protein_pool']
# find the related rxnID
idx = []
for name0 in ex_mets:
    print(name0)
    s = getRxnByReactionName(model=ecYeast, name=name0)
    if len(s) > 1:
        print("need check")
    elif len(s) == 1:
        idx.append(s[0])

# for the loop
dilutionrate = np.arange(0.05, 0.42, 0.05).tolist()
result = dict()
for k in range(len(dilutionrate)):
    print(k)
    model_tmp = ecYeast.copy()
    solution3 = DLecModelSimulate(model=ecYeast, dilution_rate=dilutionrate[k])
    fluxes = solution3.fluxes[idx]
    result[dilutionrate[k]] = list(fluxes)


# summarize the result
result_df = pd.DataFrame.from_dict(result)
result_df1 = result_df.transpose()
result_df1.columns = ex_mets
result_df1["D-glucose exchange"] = result_df1["D-glucose exchange"]*(-1)
result_df1["oxygen exchange"] = result_df1["oxygen exchange"]*(-1)
result_df2 = result_df1.drop('EX_protein_pool', axis=1) ## update

# plot
plt.figure()
sns.lineplot(x='biomass pseudoreaction', y='value', hue='variable', style="variable",
             data=pd.melt(result_df2, ['biomass pseudoreaction']))
plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0)
plt.savefig('../result/Cratree simulation based on ecModel_DLkcat.pdf', bbox_inches='tight')

