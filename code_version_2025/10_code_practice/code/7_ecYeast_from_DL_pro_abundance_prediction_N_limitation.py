# this module is mainly for ecModel simulation

from cobra.io import load_matlab_model, read_sbml_model
import sys
import matplotlib.pyplot as plt
import seaborn as sns
import pandas as pd
import numpy as np
import os


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




# compare the predicted and measured protein abundances
dir2 = "../data/ecYeast_DL_update_some_kcat.xml"
ecYeast = read_sbml_model(dir2)

# solve the model
solution3 = DLecModelSimulate(model=ecYeast, dilution_rate=0.35)

flux_max = solution3.fluxes
result = pd.DataFrame({'rxnID':flux_max.index, 'flux':flux_max.values})
result = result[result['rxnID'].str.contains("prot_")]
result['geneID'] = result['rxnID'].str.replace("prot_", "")




# input the measured protein abundances
omics_tao2 = pd.read_excel("../data/Omics_from_tao_scale.xlsx")
condition = ['gene','prot.19', 'prot.20', 'prot.21']
omics_select = omics_tao2[condition]
omics_select['average'] = omics_select.drop('gene', axis=1).apply(lambda x: x.mean(), axis=1)

# compare the predict with measured
result['pro_measured'] = singleMapping(omics_select["average"], omics_select["gene"], result['geneID'])
result = result[~result["pro_measured"].isna()]

# output the result
result.to_excel("../result/predicted_vs_measured using ecModels from DLkcat.xlsx")


# method1 absolute protein abundance mmol protein/gDW
plt.figure()
sns.regplot(x=np.log10(result['pro_measured']), y=np.log10(result['flux']), fit_reg=False)
plt.xlim(-11, 0)
plt.ylim(-11, 0)
plt.xlabel("log10(Measured_protein_level)")
plt.ylabel("log10(Predicted_protein_usage)")

from scipy.stats import pearsonr
result1 = result[result['flux'] > 0]
result1 = result1[result1['pro_measured'] > 0]
corr, ss = pearsonr(np.log10(result1['pro_measured']), np.log10(result1['flux']))
print("Correlation coefficient:", corr)
print("Correlation p_value:", ss)



# method2 calculate the RMSE
from sklearn.metrics import mean_squared_error
import math
y_actual = np.log10(result1['pro_measured'])
y_predicted = np.log10(result1['flux'])
MSE = mean_squared_error(y_actual, y_predicted)
RMSE = math.sqrt(MSE)
print("Root Mean Square Error:\n")
print(RMSE)



