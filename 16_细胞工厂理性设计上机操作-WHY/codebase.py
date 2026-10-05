# -*- coding: utf-8 -*-
import pandas as pd
import numpy as np
from cobra.flux_analysis import pfba

############################ (ec)FSEOF functions ############################
def k_matrix_filter(model, k_matrix, alpha, tol):
    # filter1.Take out rxns with no grRule:
    withGR = [rxn.id for rxn in model.reactions if rxn.gene_reaction_rule != '']
    k_matrix = k_matrix.loc[withGR,:]

    # filter2.remove out rxns that are always zero -> k=0/0=NaN:
    # all_nanList = [rxn for rxn in k_matrix.index if np.all(np.isnan(k_matrix.loc[rxn, :]))]
    all_nanList = k_matrix.apply(lambda x: np.all(pd.to_numeric(x, errors='coerce').isna()), axis=1)
    print('there are %d reactions with all NaNs' % np.sum(all_nanList))
    k_matrix = k_matrix.loc[~all_nanList]

    # 3.Replace remaining NaNs with 1s:
    # k_matrix[np.isnan(k_matrix)] = 1
    k_matrix = k_matrix.fillna(1)

    # 4.Replace any Inf value with 1000:
    k_matrix[np.abs(k_matrix) > 1000] = 1000

    # filter5.Filter out values that are inconsistent at different alphas:
    distinct_down = np.sum(k_matrix <= (1 - tol), axis=1) >= len(alpha)-2
    distinct_up = np.sum(k_matrix >= (1 + tol), axis=1) >=  len(alpha)-2
    # Identify those reactions with mixed patterns
    incons_rxns = distinct_down + distinct_up == 0
    print('there are %d reactions with mixed patterns' % np.sum(incons_rxns))
    k_matrix = k_matrix.loc[incons_rxns[~incons_rxns].index.tolist(), :]

    return k_matrix


def flux_calculation(model, targetID,c_source,c_uptake=1,model_type='GEM',tol_ratio=0.01):
    """
    Calculate the flux distribution by setting the specific objective.

    Args:
    - model: (cobra.Model) the model to optimize
    - target: (str) reaction ID for the objective reaction
    - c_source: (str) reaction ID for the main carbon source uptake reaction
    - c_uptake: (float) uptake rate for the main carbon source (default: 1)
    - tol: (float) numerical tolerance for fixing bounds(default: 1e-6)
    - prot_conc_output: (bool) whether to output protein concentration (default: False)
    - model_type: (str) type of model to optimize:'GEM'/'ecGEM' (default: 'GEM')

    Returns:
    - sol: (cobra.Solution) the simulation result

    """

    # Use pFBA to calculate the optimal flux distribution in GEMs
    if model_type=='GEM':
        model.reactions.get_by_id(c_source).bounds = -c_uptake , -c_uptake
        sol=pfba(model,objective=targetID,fraction_of_optimum=0.95)

    # Use ecFBA to calculate the optimal flux distribution in ecGEMs
    elif model_type=='ecGEM':
        model.reactions.get_by_id(c_source).bounds = c_uptake, c_uptake

        # 1.set the target objective
        model.reactions.get_by_id(targetID).bounds = 0, 1000
        model.objective = targetID
        model.objective_direction = 'max'
        sol0 = model.optimize()
        max_obj = sol0.objective_value
        # 2.Fix optimal value for objective and minimize total protein usage
        with model:
            model.reactions.get_by_id(targetID).bounds = max_obj * (1 - tol_ratio), max_obj * (1 + tol_ratio)
            prot_pool_rxnID = 'prot_pool_exchange'
            model.objective = prot_pool_rxnID
            model.objective_direction = 'min'
            sol = model.optimize()

    return sol

def flux_scanning(model, targetID, c_source,c_uptake, alpha, filterG=False,model_type='etfl',tol=0.001):
    """
    Args:
        model : ecModel/GEM model.
        targetID (str): Rxn ID for the production target reaction, a exchange reaction is recommended.
        c_source (str): Rxn ID for the main carbon source uptake reaction
        c_uptake(float): carbon source uptake rate
        alpha (np.array): scalling factor for production yield for enforced objective limits
        tol (float): numerical tolerance for fixing bounds
        filterG (bool, optional): TRUE if genes K_scores results should be filtered according to the alpha vector distribution
        model_type (str, optional): 'GEM' or 'ecGEM' model type

    Returns:
        dict: FC
            flux_WT (pd.Series): simulated flux for the wild-type (100% growth) condition
            alpha (np.array): scaling factors for production yield
            rxns (list): reaction information for each reaction in the model
            rxnNames (pd.Series): short name for each reaction in the model
            v_matrix (pd.DataFrame): matrix with simulated flux for each reaction and alpha value
            k_matrix (pd.DataFrame): matrix with normalized flux (k-score) for each reaction and alpha value
            grMatrix (pd.DataFrame): matrix with gene-reaction rules for each reaction
            k_rxns (pd.Series): median k-score across alpha values for each reaction
            genes (pd.Series): remaining genes after filtering inconsistent scores
            geneNames (pd.Series): short name for remaining genes
            k_genes (pd.DataFrame): k-score for each remaining gene
    """

    if filterG is None:
        filterG = False

    if c_uptake is None:
        c_uptake = 1

    gluc_MW=0.180156 #g/mmol
    tol_ratio=0.01

    #step 1: build reactions k_matrix
    # Simulate WT (100% max growth):
    FC = {}
    gr_rxnID = 'r_2111'

    wt_sol = flux_calculation(model=model,
                                   targetID=gr_rxnID,
                                   c_source=c_source,
                                   c_uptake=c_uptake,
                                   model_type=model_type,
                                   tol_ratio=tol_ratio)
    FC['flux_WT'] = wt_sol.fluxes
    # max_growth = FC['flux_WT'][gr_rxnID]
    # print('simulate WT-like:',max_growth)

    # simulate production in different suboptimal growth rate conditions
    FC['alpha'] = alpha
    # initialize fluxes and K_scores matrices
    rxnIDlist= [rxn.id for rxn in model.reactions]
    v_matrix = pd.DataFrame(index=rxnIDlist, columns=alpha)
    k_matrix = pd.DataFrame(index=rxnIDlist, columns=alpha)
    # simulate fluxes distribution in different growth rate conditions
    for i in range(len(alpha)):
        biomass_yield=alpha[i]
        growth=biomass_yield*c_uptake*gluc_MW
        model.reactions.get_by_id(gr_rxnID).bounds= growth*(1-tol_ratio), growth
        product_sol = flux_calculation(model=model,
                                        targetID=targetID,
                                        c_source=c_source,
                                        c_uptake= c_uptake,
                                        model_type=model_type,
                                        tol_ratio=tol_ratio)
        FC['flux_MAX'] = product_sol.fluxes
        print('when growth rate is %s,production rate is %s'%(growth,product_sol.fluxes[targetID]))
        v_matrix.iloc[:, i] = FC['flux_MAX']
        k_matrix.iloc[:, i] = FC['flux_MAX'] / FC['flux_WT']

    # step 2: calculate reactions k scores
    # filter rxns according to k_matrix
    k_matrix_filtered = k_matrix_filter(model, k_matrix, alpha, tol)
    v_matrix_filtered = v_matrix.loc[k_matrix_filtered.index.tolist(), :]
    FC["v_matrix"] = v_matrix_filtered
    FC["k_matrix"] = k_matrix_filtered
    FC["rxns"] = k_matrix_filtered.index.tolist()
    # get median k_score for target rxn
    k_rxns=pd.Series(np.mean(k_matrix_filtered, axis=1),index=k_matrix_filtered.index.tolist())
    FC["k_rxns"] = k_rxns
    # Order from highest to lowest median k_score (across alphas)
    order = np.argsort(FC["k_rxns"])[::-1]
    FC["k_rxns"] = FC["k_rxns"][order]
    FC["rxns"] = FC['k_rxns'].index.tolist()
    FC['rxnNames'] = pd.Series([model.reactions.get_by_id(rxn).name for rxn in FC['rxns']], index=FC['rxns'])
    FC["v_matrix"] = FC["v_matrix"].loc[FC["rxns"], :]
    FC["k_matrix"] = FC["k_matrix"].loc[FC["rxns"], :]

    # step 3: calculate gene k scores
    # get gene scores according k_rxns
    # build gene-rxn matrix accoreding target rxns
    grMatrix = pd.DataFrame(index=[gene.id for gene in model.genes], columns=FC['rxns'])
    for rxn in FC['rxns']:
        related_genes = [gene.id for gene in model.reactions.get_by_id(rxn).genes]
        grMatrix[rxn] = np.array([gene.id in related_genes for gene in model.genes])
    genes_rxns_dict = {}
    k_genes=pd.Series()
    for gene in grMatrix.index.tolist():
        related_rxns= grMatrix.loc[gene, grMatrix.loc[gene, :] == True].index.tolist()
        if len(related_rxns) > 0:
            related_rxns_k=FC['k_rxns'].loc[related_rxns]
            # filter genes with inconsistent k_scores
            # check if all k_scores are below 1 or all above 1
            if np.sum(related_rxns_k <= (1 - tol)) == len(related_rxns) or np.sum(related_rxns_k >= (1 + tol)) == len(related_rxns):
                genes_rxns_dict[gene] = related_rxns
                k_genes[gene] = np.mean(k_matrix_filtered.loc[related_rxns, :].values.flatten())
            else:
                print('inconsistent k_scores for gene %s'%gene)
    # Order from highest to lowest median k_score (across alphas)
    order = np.argsort(k_genes)[::-1]
    k_genes = k_genes[order]
    FC["k_genes"] = k_genes
    FC["genes"]=FC["k_genes"].index.tolist()
    FC["geneNames"] = pd.Series([model.genes.get_by_id(geneID).name for geneID in FC["genes"]],index=FC["genes"])
    grMatrix= grMatrix.loc[FC["genes"], FC["rxns"]]
    FC["grMatrix"] = grMatrix

    # step 4(optional): # Filter any value between mean(alpha) and 1
    if filterG:
        unchanged = FC['k_genes'][(FC['k_genes'] >= np.mean(alpha) - tol) & (FC['k_genes'] <= 1 + tol)].index.tolist()
        # remove unchanged genes
        FC['genes'] = [g for g in FC['genes'] if g not in unchanged]
        FC['geneNames'] = FC['geneNames'][FC['genes']]
        FC['k_genes'] = FC['k_genes'][~unchanged]
        FC['grMatrix'] = FC['grMatrix'].loc[FC['genes'], :]
        # Update results for gene-related reactions (remove remaining reactions without any associated gene in rxnGeneM)
        FC['rxns'] = FC['grMatrix'].columns[FC['grMatrix'].any(axis=0)].tolist()
        FC['v_matrix'] = FC['v_matrix'].loc[FC['rxns'], :]
        FC['k_matrix'] = FC['k_matrix'].loc[FC['rxns'], :]
        FC['k_rxns'] = FC['k_rxns'][FC['rxns']]
        FC['rxnNames'] = [model.reactions.get_by_id(rxnID) for rxnID in FC['rxns']]

    return FC


def run_FSEOF(model, targetID, c_source,c_uptake, alphaLims, Nsteps,model_type='GEM'):
    """
    Run Flux-scanning with Enforced Objective Function for a specified production target.

    Arguments:
    - model:
        GEM/ecGEM model
    - targetID: str
        Reaction ID for the production target reaction, a exchange reaction is recommended
    - c_source: str
        Reaction ID for the main carbon source uptake reaction (make sure that the correct directionality is indicated)
    - alphaLims: tuple of float
        Minimum and maximum biomass yield [gDw/mmol Csource] for enforced objective limits
    - Nsteps: int
        Number of steps for suboptimal objective in FSEOF
    - model_type: str
        Type of model (GEM or ecGEM)
    Returns:
    - dict
        Dictionary with the results of the FSEOF analysis
    """

    # Define alpha vector for suboptimal enforced objective values
    alphaV = np.linspace(alphaLims[0], alphaLims[1], Nsteps)
    # Run FSEOF analysis
    results = flux_scanning(model=model, targetID=targetID, c_source=c_source,c_uptake=c_uptake, alpha=alphaV,model_type=model_type)

    # Create gene table
    geneTable = pd.DataFrame(index=results['genes'], columns=["gene_names", "k_score"], dtype=float)
    geneTable.index.name = "gene_IDs"
    geneTable["gene_names"] = results["geneNames"]
    geneTable["k_score"] = results['k_genes']
    results["geneTable"] = geneTable

    # Create reaction table
    rxnTable = pd.DataFrame(index=results['rxns'], columns=["rxn_names", "k_score"], dtype=float)
    rxnTable.index.name = "rxn_IDs"
    rxnTable["rxn_names"] = results["rxnNames"]
    rxnTable["k_score"] = results['k_rxns']
    results["rxnTable"] = rxnTable

    # Remove redundant output fields
    results.pop('k_rxns', None)
    results.pop('rxns', None)
    results.pop('rxnNames', None)
    results.pop('genes', None)
    results.pop('geneNames', None)
    results.pop('k_genes', None)
    results.pop('flux_MAX', None)

    return results


############################ optForce functions ############################
def run_FVA(model, target_rxn, target_const_rate, flux_constraints={}):
    with model as m:
        m.objective = target_rxn
        m.objective_direction = 'max'
        for rxn_id in flux_constraints:
            m.reactions.get_by_id(rxn_id).bounds = flux_constraints[rxn_id]
        target_max = m.slim_optimize()

    fva_result = {}
    for rxn in model.reactions:
        with model as m:
            tmp = target_const_rate * target_max
            m.reactions.get_by_id(target_rxn).bounds = (tmp, tmp)

            for rxn_id in flux_constraints:
                m.reactions.get_by_id(rxn_id).bounds = flux_constraints[rxn_id]

            m.objective = rxn.id

            m.objective_direction = 'min'
            min_flux = m.slim_optimize()
            m.objective_direction = 'max'
            max_flux = m.slim_optimize()
            fva_result[rxn.id] = (min_flux, max_flux)
    return fva_result


def compare_FV_range(bio_result, target_result):
    down_targets = {}
    up_targets = {}
    for rxn_id in target_result:
        target_lb, target_ub = target_result[rxn_id]
        bio_lb, bio_ub = bio_result[rxn_id]
        if target_lb > bio_ub:
            up_targets[rxn_id] = {'bio': (bio_lb, bio_ub),
                                  'target': (target_lb, target_ub)}
        elif target_ub < bio_lb:
            down_targets[rxn_id] = {'bio': (bio_lb, bio_ub),
                                    'target': (target_lb, target_ub)}
    return down_targets, up_targets


def recordResults(target_rxn, down_targets, up_targets):
    if not os.path.exists('./output'):  # If the directory does not exist, make it
        os.makedirs('./output')

    with open(f'./output/OptForce_{target_rxn}.txt', 'w') as fp:
        fp.write(f'Target reaction: {target_rxn}\n')
        fp.write('\nDown targets\tBio lb\tBio ub\tTarget lb\tTarget ub\n')
        for rxn_id in down_targets:
            bio_lb, bio_ub = down_targets[rxn_id]['bio']
            target_lb, target_ub = down_targets[rxn_id]['target']
            fp.write(f'{rxn_id}\t{bio_lb:4f}\t{bio_ub:4f}\t{target_lb:4f}\t{target_ub:4f}\n')
        fp.write('\nUp targets\tBio lb\tBio ub\tTarget lb\tTarget ub\n')
        for rxn_id in up_targets:
            bio_lb, bio_ub = up_targets[rxn_id]['bio']
            target_lb, target_ub = up_targets[rxn_id]['target']
            fp.write(f'{rxn_id}\t{bio_lb:4f}\t{bio_ub:4f}\t{target_lb:4f}\t{target_ub:4f}\n')


def rxns_to_genes(model, rxnIDlist):
    geneIDlist = list()
    for rxnID in rxnIDlist:
        for g in model.reactions.get_by_id(rxnID).genes:
            geneIDlist.append(g.id)

    geneIDlist = list(set(geneIDlist))
    return geneIDlist


def parse_optforce_result(model,down_rxns,up_rxns):
    # transform reaction targets to gene targets
    down_targets_idlist = down_rxns.keys()
    up_targets_idlist = up_rxns.keys()
    down_gene_idlist = rxns_to_genes(model, down_targets_idlist)
    up_gene_idlist = rxns_to_genes(model, up_targets_idlist)

    # remove conflicting targets
    common_gene_idlist = list(set(down_gene_idlist) & set(up_gene_idlist))
    if len(common_gene_idlist) > 0:
        down_gene_idlist = list(set(down_gene_idlist) - set(common_gene_idlist))
        up_gene_idlist = list(set(up_gene_idlist) - set(common_gene_idlist))

    df_gene_result = pd.DataFrame(down_gene_idlist + up_gene_idlist, columns=['geneID'])
    df_gene_result.loc[df_gene_result['geneID'].isin(down_gene_idlist), 'action'] = 'KD'
    df_gene_result.loc[df_gene_result['geneID'].isin(up_gene_idlist), 'action'] = 'OE'

    return df_gene_result