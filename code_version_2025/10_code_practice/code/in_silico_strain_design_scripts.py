# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

# this script will be not updated. Refer to code from repo - strain design



# Import packages
import math
from cobra.manipulation import *
import pandas as pd

# function for the in silico strain design
def find(rxn_name_list, specific_string, equal=False):
    '''
    function to find the index of element in a list contains a specific string, refer the function
    from matlab
    :param rxn_name_list:
    :param specific_string:
    :param equal: if true will return the indexes in rxn_name_list which is equal to specific_string
    :return: index of the element from the list contains the specific string
    example:
    rxn_name_list = ['a_b','e_b','c','d']
    specific_string = '_b'
    '''
    s0 = rxn_name_list
    if equal ==False:
        index =[i for i,x in enumerate(s0) if specific_string in x]
    else:
        index = [i for i, x in enumerate(s0) if specific_string == x]
    return index


# function
def ecYeastMinimalMedia(model):
    """
    This function is used to define a simple media for ecYeast
    :param model:
    :return: a model with the defined the minimal media
    """
    rxnID = []
    rxnName = []
    for i, x in enumerate(model.reactions):
        rxnID.append(x.id)
        rxnName.append(x.name)

    exchange_rxn =[x for x, y in zip(rxnID, rxnName) if '_REV' in x and 'exchange' in y]
    # first block any uptake
    for i, x in enumerate(exchange_rxn):
        rxn0 = exchange_rxn[i]
        print(rxn0)
        model.reactions.get_by_id(rxn0).upper_bound = 0

    #Allow uptake of essential components
    model.reactions.get_by_id("r_1654_REV").upper_bound = 10000 #ammonium exchange (reversible)
    model.reactions.get_by_id("r_1861_REV").upper_bound = 10000 #iron(2+) exchange (reversible)
    model.reactions.get_by_id("r_2100_REV").upper_bound = 10000 #water exchange (reversible)
    model.reactions.get_by_id("r_1992_REV").upper_bound = 10000 #oxygen exchange (reversible)
    model.reactions.get_by_id("r_2005_REV").upper_bound = 10000 #phosphate exchange (reversible)
    model.reactions.get_by_id("r_2060_REV").upper_bound = 10000 #sulphate exchange (reversible)
    model.reactions.get_by_id("r_1832_REV").upper_bound = 10000 #H+ exchange (reversible)
    return model


def frange(start, stop, step):
    """
    This function is like range, step can be float value, like 0.1
    :param start:
    :param stop:
    :param step:
    :return:
    """
    i = start
    while i < stop:
        yield i
        i += step
def changeInfTo1000(model):
    """
    This function is used to change the inf in the uppper bound of reaction into 10000
    :param model:
    :return:
    """
    # update the upper bound
    for x in model.reactions:
        print(x.id)
        if x.upper_bound == math.inf:
            model.reactions.get_by_id(x.id).upper_bound = 1000
        else:
            pass
    return model

# This function obtain the mutant model after removing one gene
def getModelWithRemoveGene (model0, gene_remove0):
    '''

    :param model0:
    :param gene_remove0: should be a list
    :return:
    '''
    model1 = model0.copy()
    remove_genes(model1, gene_remove0, remove_reactions=True) # this script can't be used
    return model1


# function to carry out the sensitivity analysis
# here we firstly increase all the upper bound and test each one
def sensitivityEnzymeLevel(Model, rxn_list, rxnID, down_factor = 1, up_factor = 10):
    with Model:
        # exclude the rxnID from rxn_list
        rxn_list_test = [i for i in rxn_list if i != rxnID]
        # open all other protein usage
        for i in rxn_list_test:
            model.reactions.get_by_id(i).bounds = (0, 10)
        # set the new upper bound for the test reaction
        upper_bound = Model.reactions.get_by_id(rxnID).upper_bound
        # decrease and increase the upper bound
        upper_down = upper_bound * down_factor
        upper_up = upper_bound * up_factor

        # calculate the growth when adjust the bounds
        Model.reactions.get_by_id(rxnID).upper_bound = upper_down
        solution1 = Model.optimize()
        print('solution for the down-regulated:')
        Model.summary()
        growth1 = solution1.objective_value

        Model.reactions.get_by_id(rxnID).upper_bound = upper_up
        solution2 = Model.optimize()
        print('solution for the up-regulated:')
        Model.summary()
        growth2 = solution2.objective_value

    fold_change = growth2 / growth1
    print('growth increased by:')
    print(fold_change)
    return growth1, growth2, fold_change




# here we use a function to get the proper upper bound under the fixed growth rate
# the function was used to get a proper upper bound to make the cell have a normal growth rate
def FlexProtein(Model,rxn_list, constraint, biomassID='r_2111'):
    with Model:
        # exclude the rxnID from rxn_list
        # open all other protein usage
        for i in rxn_list:
            Model.reactions.get_by_id(i).bounds = (0, 10)

        # set the bounds for specific reactions with measured values
        for key in constraint:
            Model.reactions.get_by_id(key).bounds = (constraint[key], constraint[key])

        solution= Model.optimize()
        Model.summary()
        # fix the growth rate
        Model.reactions.get_by_id(biomassID).bounds = (0.99*solution.objective_value, 1.01*solution.objective_value)

        # minimization of all the proteins usage as a whole
        for rxn in Model.reactions:
            print(rxn.id)
            rxn.objective_coefficient = 0
        print(Model.objective._expression)
        for rxn in rxn_list:
            print(rxn)
            Model.reactions.get_by_id(rxn).objective_coefficient = -1
        print(Model.objective._expression)
        sol3 = Model.optimize()
        #Model.summary()
        sol3.fluxes['r_2111']
        sol3.fluxes['r_1714_REV']
        # then get the bounds of the limited protein level
    return sol3.fluxes


# this function is used to find the prot_exchanage reactiona and its upper bound
def getProtExchangeInf(ecModel):
    rxn0 = []
    bound0 = []
    with ecModel:
        for rxn in ecModel.reactions:
            s = rxn.id
            print(s)
            if 'prot_' in s and '_exchange' in s:
                rxn0.append(s)
                bound0.append(rxn.upper_bound)
            else:
                pass
    return rxn0, bound0

# This function is used to find the limited proteins in the proteomics constrained model
# the normal_miu is the largest growth rate during the process to find the limited proteins
def findLimitedProtein(ecModel, constraint0, normal_miu):
    Model = ecModel.copy()
    new_objective = Model.problem.Objective(
            (6 * Model.reactions.r_1714_REV.flux_expression - 6 * constraint0['r_1714_REV']) ** 2 +
            (2 * Model.reactions.r_1761_REV.flux_expression - 2 * constraint0['r_1761_REV']) ** 2 +
            (Model.reactions.r_1672.flux_expression - constraint0['r_1672']) ** 2 +
            (Model.reactions.r_1992_REV.flux_expression - constraint0['r_1992_REV']) ** 2 +
            (Model.reactions.r_2111.flux_expression * 1000 / 26.9 - constraint0['r_2111'] * 1000 / 26.9) ** 2 +
            (40 * Model.reactions.DM_phytoene.flux_expression - 40 * constraint0['DM_phytoene']) ** 2 +
            (40 * Model.reactions.DM_lycopene.flux_expression - 40 * constraint0['DM_lycopene']) ** 2 +
            (40 * Model.reactions.DM_b_carotene.flux_expression - 40 * constraint0['DM_b_carotene']) ** 2,
            direction='min')

    Model.objective = new_objective
    #solution2 =Model.optimize(objective_sense=None)
    newGR = 0
    limit_protein = []
    while newGR <= 0.9999 * normal_miu:
        solution2 = Model.optimize(objective_sense=None)
        # loops to find the limited protein and then change its upper bound
        reduce_cost = solution2.reduced_costs
        reduce_cost0 = reduce_cost.sort_values(ascending=True)
        # extract the pro_exchange
        reduce_cost1 = reduce_cost0[
            [i for i, x in enumerate(reduce_cost0.index) if 'prot_' in x and 'exchange' in x]]
        limit_protein_index = list(reduce_cost1[reduce_cost1 == min(reduce_cost1)].index)[0]
        Model.reactions.get_by_id(limit_protein_index).upper_bound = 1000
        try:
            solution2 = Model.optimize(objective_sense=None)
            newGR = solution2.fluxes['r_2111']
            print(limit_protein_index, 'new growth:', newGR, 'reduce cost:', min(reduce_cost1))
            limit_protein.append(limit_protein_index)
        except:
            pass
    return limit_protein


# Find the new bound of limited protein by minimizing the sum of (new_protein_bound – measured_protein_bound)^2
# The input ecModel should be firstly constrainted by the measured rate.
def FindNewBound(ecModel,constraint, limit_proteins):
    with ecModel:
        test_rxn, upper_bound = getProtExchangeInf(ecModel=ecModel)
        ecModel.reactions.get_by_id("r_1714_REV").bounds = (0, 1000)  # open the glucose uptake rate
        for i, x in enumerate(test_rxn):
            print(x)
            ecModel.reactions.get_by_id(x).upper_bound = 10
        # here the constraint is two strict
        # in this case, the obtain the result using minimization the measured and predicted proteins will be similar to previous method
        for key in constraint:
            if key in ['r_2111', 'r_1992_REV', 'r_1672']:
                ecModel.reactions.get_by_id(key).bounds = (0.99 * constraint[key], 1.01 * constraint[key])
            else:
                ecModel.reactions.get_by_id(key).bounds = (constraint[key], constraint[key])
        solution2 = ecModel.optimize(objective_sense=None)
        all_flux = solution2.fluxes
        # find the new upper bound
        # import pandas as pd
        df0 = pd.DataFrame({'growth_l': [None] * len(test_rxn), 'growth_u': [None] * len(test_rxn),
                           'fold_change': [None] * len(test_rxn)});
        df0['protein'] = test_rxn
        df0['protein_level'] = upper_bound
        df0['new_upper_bound'] = [None] * len(test_rxn)
        new_protein = []  # note: don't update value in dataframe during the loop
        for i, x in df0.iterrows():
            print(i, x['protein'])
            new_bound = all_flux[x['protein']]
            if x['protein'] in limit_proteins:
                if df0['protein_level'][i] < new_bound:
                    new_protein.append(new_bound)
                else:
                    new_protein.append('proteomics data can be used')
            else: # in this step we need bring in more limited proteins, but it is essential for the normal prediction in the normal strain
                if df0['protein_level'][i] < new_bound:
                    new_protein.append(new_bound)
                else:
                    new_protein.append('proteomics data can be used')
        df0['new_upper_bound'] = new_protein
    return df0




# solve the proteomics constrainted model using general procedures
def SolveProteinConstraint(Model, rxn_list, constraint, biomassID ='r_2111'):
    with Model:
        for key in constraint: # for the strains with product, we need add the constraint here
            Model.reactions.get_by_id(key).bounds = (constraint[key], constraint[key])
        Model.objective = biomassID
        sol = Model.optimize()
        growth = sol.fluxes[biomassID]
        # fix the biomass
        Model.reactions.get_by_id(biomassID).bounds =(0.99*growth, 1.01*growth)

        # fix the glucose
        new_objective = Model.problem.Objective(
            180 * Model.reactions.r_1714_REV.flux_expression + 46 *
            Model.reactions.r_1761_REV.flux_expression,
            direction='min')
        Model.objective = new_objective
        sol2 = Model.optimize()
        Model.summary()
        print('substrate uptake rate:')
        print(sol2.fluxes['r_1714_REV'])
        print(sol2.fluxes['r_1761_REV'])

        # fix the GUR
        Model.reactions.get_by_id('r_1714_REV').upper_bound = 1.02*sol2.fluxes['r_1714_REV'] # using factor 1.02 to avoid the too fixed upper bound
        Model.reactions.get_by_id('r_1761_REV').upper_bound = sol2.fluxes['r_1761_REV']

        # minimization of all the proteins usage as a whole
        for rxn in Model.reactions:
            print(rxn.id)
            rxn.objective_coefficient = 0
        print(Model.objective._expression)
        for rxn in rxn_list:
            print(rxn)
            Model.reactions.get_by_id(rxn).objective_coefficient = -1
        print(Model.objective._expression)
        sol3 = Model.optimize()
        Model.summary()
        sol3.fluxes['r_1714_REV']
        sol3.fluxes['r_1761_REV']
        # then get the bounds of the limited protein level
    return sol3.fluxes














