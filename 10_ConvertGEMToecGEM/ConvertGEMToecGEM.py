#!/usr/bin/env python
# coding: utf-8

# In[1]:


import cobra
from cobra.io import load_matlab_model, read_sbml_model, load_json_model
from cobra import Model, Reaction, Metabolite
import pandas as pd
import numpy as np


# **导入.json格式模型**

# In[2]:


model_file = 'data/Halo_GEM_v1.json'
model = cobra.io.load_json_model(model_file)


# # 定义cenvertToIrrev函数

# In[3]:


def convertToIrrev(model,rxns=None):
        if rxns is None:
            rxns = model.reactions
        converted_reaction_list = []
        for rxn in rxns:
            rxn = rxn.copy()
            if rxn.lower_bound>=0:
                converted_reaction_list.append(rxn)
            else:
                rxn_REV= rxn.copy()
                rxn_REV.id= rxn.id + '_REV'
                rxn_REV.name = rxn.name + ' (reversible)'
                rxn_REV.add_metabolites({met:-rxn.get_coefficient(met)*2 for met in rxn.metabolites})
                rxn_REV.lower_bound = 0
                rxn_REV.upper_bound = -rxn.lower_bound
                rxn.lower_bound = 0
                converted_reaction_list.extend([rxn,rxn_REV])
        irrevModel = Model()
        irrevModel.add_reactions(converted_reaction_list)
        return irrevModel


# **`def convertToIrrev(model,rxns=None):`**  
#     # 定义一个函数，convert to irrevocable，将模型中的反应转化为不可逆反应。它接受两个参数：model（模型对象）和 rxns（可选的反应列表）。如果 rxns 没有被提供，它将默认为 None。在 Python 中，当你定义一个函数并为其参数指定了默认值时，如果在调用函数时没有提供该参数的值，那么函数会使用默认值。但如果在调用函数时提供了一个值，那么这个提供的值会覆盖默认值。  
#   **`    if rxns is None:`**#rxn，即reaction；rxns是reaxtions。  
#         **`        rxns = model.reactions`**#，如果rxns为None（输入的时候，未指定rxns的值），将 rxns 设置为模型中的所有反应。  
#     **`    converted_reaction_list = []`**#创建一个列表，接受转化后的反应  
#    **`    for rxn in rxns:`**  
#         **`        rxn = rxn.copy()`**#创建一个反应的副本，避免修改原模型。  
#         

#  **#irreversible reactions对于反应中的不可逆反应，执行以下操作：**  
#         `    if rxn.lower_bound>=0:`    
#             `        converted_reaction_list.append(rxn)`
#             
# #检查反应的下界是否大于或等于 0，这表示反应是不可逆的。在代谢网络模型中，每个反应都可以有一个下界（lower bound）和一个上界（upper bound），这些界限定义了反应可以进行的方向和速率。下界和上界通常用来表示反应速率的可能范围，其中：
#         
# 下界（lower_bound）表示反应可以以多快的速率向负方向（逆反应方向）进行。
#                               
# 上界（upper_bound）表示反应可以以多快的速率向正方向（正反应方向）进行。
#                               
# 在许多代谢模拟软件和库中，如 COBRApy（一个流行的代谢网络分析库），反应的下界和上界用来定义反应的可逆性：
#         
# 如果一个反应的下界大于或等于 0（lower_bound >= 0），这意味着反应不能自然地向负方向进行，即它不能自发地逆转。在这种情况下，我们说反应是“不可逆的”。
#          
# 如果一个反应的下界小于 0（lower_bound < 0），这意味着反应可以向负方向进行，即它可以自发地逆转。在这种情况下，我们说反应是“可逆的”。
#         
# 因此，当我们说“检查反应的下界是否大于或等于 0，这表示反应是不可逆的”，我们是在检查反应是否被设定为只能向正方向进行，而不能向负方向进行。这是代谢网络分析中的一个重要概念，因为不同的分析方法可能要求所有反应都是不可逆的。
#         
#  **# reversible_reacions对于反应中的可逆反应，执行以下操作：**  
#         
# **`    else:`**#如果反应的下界小于0，即反应是可逆反应。  
#             **`        rxn_REV= rxn.copy()`**#创建另一个反应副本，用于创建反向反应。  
#             **`        rxn_REV.id= rxn.id + '_REV'`**  
#             `        rxn_REV.name = rxn.name + ' (reversible)'`#设置反向反应的id和名称  
#             `        rxn_REV.add_metabolites({met:-rxn.get_coefficient(met)*2 for met in rxn.metabolites})`#更新反向反应中代谢物的系数，使其与原反应相反，并乘以 2 以保持平衡。  
#             `        rxn_REV.lower_bound = 0`#设置逆反应的下界为0  
#             `        rxn_REV.upper_bound = -rxn.lower_bound`#设置逆反应的上界为原反应下界的相反数（相当于把原本的反应拆分为一正一负两个反应）  
#             `        rxn.lower_bound = 0`#将正向反应的下界修改为0，变为不可逆反应  
#             `        converted_reaction_list.extend([rxn,rxn_REV])`#将原反应和反向反应都添加到 converted_reaction_list 列表中。  

# **build irrevModel完成可逆反应向不可逆反应的转换后，创建一个新的不可逆模型来接收这些反应。**

# `    irrevModel = Model()  `
# 
# #创建一个新的模型类，cobra模型类中通常包含反应（Reactions）、代谢物（Metabolites）、反应的边界（Reaction Bounds）、反应的基因关联（Gene Associations）、目标函数（Objective Function）：模型可以定义一个或多个目标函数，用于优化分析。目标函数通常是一个线性组合的代谢物通量，用于最大化或最小化特定的生物学目标、约束（Constraints）：模型包含一组约束，确保反应的进行遵守质量守恒和热力学定律。  
# 模型参数（Model Parameters）：模型可能包含其他参数，如生长速率、通量控制系数等，用于模拟和分析、解决方案（Solution）：在进行优化后，模型会包含一个解决方案对象，其中包含优化结果，如最大生长速率、反应通量分布等。模型 ID 和名称（Model ID and Name）和其他属性（Other Attributes）。  
#     `    irrevModel.add_reactions(converted_reaction_list)`  
#     `    return irrevModel`#convertToIrrev函数的返回值为一个不可逆反应的Model类

# # 转为不可逆模型

# In[4]:


irrModel=convertToIrrev(model)


# 设置上下限

# In[5]:


def set_bound(model,rxn_id,lb=None,ub=None,eq=None):
    rxn = model.reactions.get_by_id(rxn_id)
    if lb is not None: rxn.lower_bound = lb
    if ub is not None: rxn.upper_bound = ub
    if eq is not None: 
        rxn.lower_bound = eq
        rxn.upper_bound = eq
set_bound(irrModel,'Exchange_Glucopyranose',lb=0,ub=10)


# ## 打印模型查看模型反应数变化

# In[6]:


def report_model(model):
    print('Reactions:',len(model.reactions))
    print('Metabolites:',len(model.metabolites))


# In[7]:


report_model(irrModel)


# In[8]:


report_model(model)


# ## 验证模型对生物量和产量的预测准确性

# In[9]:


def set_yeast_extraction(model,ub=1000,lb=-1000):
    amino_acids = ['CYS_c','L-ASPARTATE_c','GLT_c','PHE_c','GLY_c',
                  'HIS_c','ILE_c','LYS_c','LEU_c','MET_c','ASN_c',
                  'PRO_c','GLN_c','ARG_c','SER_c','THR_c','VAL_c',
                  'TRP_c','L-ALPHA-ALANINE_c','TYR_c','L-LACTATE','CPD-15990', 'ADENINE']
    for met_id in amino_acids:
        exg_rxn = model.reactions.get_by_id('Exchange_'+met_id.replace('_c',''))
        exg_rxn.upper_bound = ub
        exg_rxn.lower_bound = lb


# #设置酵母与外界交换氨基酸的反应的上下界（-1000~1000）。

# In[10]:


def test_biomass_production(model,show=False):
    with model:
        set_yeast_extraction(model,lb=0,ub=0)
        model.objective = 'Biomass_v1'
        model.objective_direction = 'max'
        s1 = model.optimize()
        if show:print('Without YE:',s1)

    with model:
        set_yeast_extraction(model,lb=0,ub=1000)
        model.objective = 'Biomass_v1'
        model.objective_direction = 'max'
        s2 = model.optimize()
        if show:print('With unlimited YE:',s2)
    
    if not show:return s1,s2


# #模拟模型的产Biomass反应，返回两个值：第一个为完全阻断氨基酸摄取时(Without YE)的Biomass产量，第二个为允许氨基酸摄取时(With unlimited YE)的Biomass产量。

# In[11]:


test_biomass_production(irrModel,True)


# # 添加kcat约束

# In[12]:


import pandas as pd


# ## 导入kcat参数

# In[13]:


df_enz_kcat = pd.read_csv('data/mapped_kcats.csv',index_col=0)


# #读取kcat文件为df_enz_kcat这个DataFrame。

# In[14]:


def prepare_kcats_dict(irrModel, df_enz_kcat,col='log10_kcat_mean'):
    kcat_dict = {}
    for rxn_id in df_enz_kcat.index:
        rxn = irrModel.reactions.get_by_id(rxn_id)
        for gene in rxn.genes:
            kcat_dict[(gene.id,rxn_id)] = 10**df_enz_kcat.loc[rxn_id,col]*3600 # kcat in 1/h
    return kcat_dict


# #定义prepare_kcats_dict函数，从df_enz_kcat中提取与基因和反应id一一对应的kcat值为字典kcat_dict，格式为：
# 
# (enz_id,rxn_id)：kcat 
# 
# 例如：
# 
#     kcat_dict = {
#                      ('E1','R3'): 10,
#                      ('E1','R4'): 100,
#                      ('E2','R4'): 0.10,
#                      ('E3','R5'): 5,
#                      ('E3','R5_REV'): 12,
#                      ('E4','R5'): 15,
#                      ('E4','R5_REV'): 90,
#                   }
#                   
# 此外，kcat值文件mapped_kcats.csv中，kcat值以log10_kcat格式存储，所以prepare_kcats_dict函数还将该值转化为单位为1/h的kcat值，即
# 
# kcat=10^log10_kcat_mean*3600        

# In[15]:


kcat_dict = prepare_kcats_dict(irrModel,df_enz_kcat,'log10_kcat_max')


# ## 定义添加kcat参数的函数

# In[16]:


def addEnzymesToRxn(rxn, kcat, protIDs, rxn_index=None):
    
    e_rxn      = rxn.copy()
    if rxn_index is not None:
        e_rxn.id   = e_rxn.id + 'No{0}'.format(rxn_index)
        e_rxn.name = e_rxn.name + ' (No{0})'.format(rxn_index)
    prots = [item.strip() for item in protIDs.split('and')]
    
    
    # get compartment
    comp = None
    for met in rxn.metabolites:
        comp = met.compartment
        if rxn.get_coefficient(met)<0: comp = met.compartment
    
    # create Metabolite object for each protein and create exchange reaction
    prot_mets = []
    prot_exchange_rxns = []
    for prot in prots:
        prot_met = Metabolite('prot_{0}_{1}'.format(prot,comp))
        prot_met.compartment =  comp
        prot_mets.append(prot_met)
        
        # add excange reaction of protein
        excg_rxn = Reaction('prot_{0}_exchange'.format(prot))
        excg_rxn.lower_bound = 0
        excg_rxn.gene_reaction_rule = prot
        excg_rxn.add_metabolites({prot_met:1})
        prot_exchange_rxns.append(excg_rxn)
        
    # add enzymes into reaction
    e_rxn.add_metabolites({prot_met:-1./kcat for prot_met in prot_mets})
    e_rxn.gene_reaction_rule = protIDs

    return e_rxn, prot_exchange_rxns


# #将酶作为代谢物添加到特定的反应中，并为每个酶创建一个交换反应。这样，可以在代谢模型中考虑酶的动力学特性，例如通过 kcat 值来调节反应的速率。此外，函数还允许为同工酶（isoenzymes）创建不同的反应版本，通过 rxn_index 参数来区分。

# In[17]:


def parse_gr_rule(gr):
    complexes = [item.strip().replace('(','').replace(')','') for item in gr.split('or')]
    if len(complexes) < 2 and len(complexes[0]) < 1:
        complexes = []
    return complexes


# #将基因反应规则中的复杂表达式分解成单独的组件，便于后续处理基因反应规则。
# 
#  For example: 
#     
#     Input         : Output
#     A or B        : ["A", "B"]
#     (A and B)     : ["A and B"]
#     (A and B) or C: ["A and B","C"]

# In[18]:


def getArmReaction(rxn):
    # 1. create intermediate metabilite
    rxnID = rxn.id
    comp = None
    for met in rxn.metabolites:
        comp = met.compartment
        if rxn.get_coefficient(met)<0: comp = met.compartment
    
    pmet = Metabolite('pmet_{0}'.format(rxnID),compartment=comp)
    
    # 2. create arm reaction: 
    arm_rxn                    = Reaction('arm_{0}'.format(rxnID))
    arm_rxn.name               = rxn.name + ' (arm)'
    arm_rxn.subsystem          = rxn.subsystem
    arm_rxn.lower_bound        = rxn.lower_bound  
    arm_rxn.upper_bound        = rxn.upper_bound  
    arm_rxn.gene_reaction_rule = ''
    
    mets = {met:rxn.get_coefficient(met) for met in rxn.metabolites if rxn.get_coefficient(met)<0}
    mets[pmet] = 1
    
    arm_rxn.add_metabolites(mets)
    
    # 3. change orignal reaction to pmet --> C + D 
    rxn_new = rxn.copy()
    rxn_new.subtract_metabolites({met:rxn_new.get_coefficient(met) for met in rxn_new.metabolites if rxn_new.get_coefficient(met)<0})
    rxn_new.add_metabolites({pmet:-1})
    
    return rxn_new, arm_rxn


# #将一个代谢反应转换为一个两步过程，其中原始反应的底物被一个中间代谢物（pmet）替代，而原始产物则保留在修改后的反应中。
# 
#     Original reaction: A + B --> C + D
#     Arm reaction    : A + B --> pmet   (no gr rule)
#     Change the orginial reaction to:  pmet --> C + D (use old gr rule)

# In[19]:


def convertToEnzymeModel(model,kcats):
    converted_reaction_list = []
    protein_exchange_rxns = {}
    for rxn in model.reactions:
        complexes = parse_gr_rule(rxn.gene_reaction_rule)
        
        # 1. for those reactions without genes 
        if len(complexes) <1: 
            converted_reaction_list.append(rxn)
            continue
        
        # 2. for those reactions with genes, but no kcat
        first_gene = [gene.id for gene in rxn.genes][0]
        if kcats.get((first_gene,rxn.id),None) is None:
            converted_reaction_list.append(rxn)
            continue
        
        # 3. for those reactions with isoenzymes, add arm reaction
        if len(complexes) >1:
            rxn_new, arm_rxn = getArmReaction(rxn)
            converted_reaction_list.append(arm_rxn)
            
            for i,complx in enumerate(complexes):
                prots = [item.strip() for item in complx.split('and')]
                kcat = kcats[(prots[0],rxn.id)]
                e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn_new, kcat, complx,rxn_index=i+1)
                
                converted_reaction_list.append(e_rxn)
                for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rxn.id] = prot_exchange_rxn
                
            continue
         
        if len(complexes) == 1:
            complx = complexes[0]
            prots = [item.strip() for item in complx.split('and')]
            kcat = kcats[(prots[0],rxn.id)]
            e_rxn, prot_exchange_rxns = addEnzymesToRxn(rxn, kcat, complx,rxn_index=1)
            converted_reaction_list.append(e_rxn)
            for prot_exchange_rxn in prot_exchange_rxns: protein_exchange_rxns[prot_exchange_rxn.id] = prot_exchange_rxn
    
    eModel = Model()
    eModel.add_reactions(converted_reaction_list)
    eModel.add_reactions(protein_exchange_rxns.values())
    eModel.enzymes = set([exgrxn.split('_')[1] for exgrxn in protein_exchange_rxns.keys()])
    print('Number of enzymes:', len(eModel.enzymes))
    
    return eModel


# #定义了一个名为 convertToEnzymeModel 的函数，它用于将一个基础的代谢模型转换为包含酶动力学参数的模型。
# 函数考虑三种情况：没有基因的反应、有基因但没有 kcat 值的反应和有同工酶的反应。其中，对于有同工酶的反应，创建ArmReaction。
# 

# In[20]:


eModel = convertToEnzymeModel(irrModel,kcat_dict)


# ### 打印模型查看模型反应数变化

# In[21]:


report_model(eModel)


# ### 验证模型对生物量和产量的预测准确性
# 

# In[22]:


test_biomass_production(eModel)


# # 绝对定量蛋白质组学数据添加

# In[23]:


dfomics = pd.read_csv('data/protein_abundance.csv',index_col=0)
dftot = pd.read_csv('data/total_protein_abundance.csv',index_col=0)
dfmws = pd.read_csv('data/protein_mws.csv',index_col=0)


# In[24]:


dfmws = dfmws/1000
dfmws.head()


# In[25]:


def prepare_omics_for_one_condition(dfomics,dftot,dfmws,condition_id,enzymes):

    non_measured = []
    measured = {}
    
    assert condition_id in dfomics.columns
    for enz in enzymes:
        try: ab = dfomics.loc[enz,condition_id]
        except: ab = np.nan
        
        if np.isnan(ab): non_measured.append(enz)
        else: measured[enz] = ab
    
    tot_measured_proteins = np.sum([0 if np.isnan(dfomics.loc[ind,condition_id]) 
                                    else dfomics.loc[ind,condition_id]*dfmws.loc[ind,'MW'] 
                                    for ind in dfomics.index])
    
    tot_measured_enzymes  = np.sum([measured[ind]*dfmws.loc[ind,'MW'] for ind in measured.keys()])
    tot_non_measured      = dftot.loc[condition_id,'Ptot'] - tot_measured_proteins
    enzyme_fraction       = tot_measured_enzymes/tot_measured_proteins
    prot_pool = tot_non_measured*enzyme_fraction
    
    print('tot_proteins         :',dftot.loc[condition_id,'Ptot'])
    print('tot_measured_proteins:',tot_measured_proteins)
    print('tot_measured_enzymes :',tot_measured_enzymes)
    print('tot_non_measured     :',tot_non_measured)

    print('enzyme_fraction      :',enzyme_fraction)
    print('prot_pool            :',prot_pool)
    
    return measured, non_measured, prot_pool,enzyme_fraction


# #确定哪些酶的丰度是已知的，哪些是未知的，并计算未测量酶的总丰度

# In[26]:


def prepare_omics_for_one_condition(dfomics, dftot, dfmws, condition_id, model):
    non_measured = []  # 未测量的酶
    measured = {}      # 已测量的酶
    
    # 检查condition_id是否在dfomics中
    assert condition_id in dfomics.columns

    # 遍历模型中的酶
    for enz in model.enzymes:
        try:
            ab = dfomics.loc[enz, condition_id]  # 获取该酶在条件下的表达量
        except KeyError:  # 如果没有找到该酶，赋值为NaN
            ab = np.nan
        
        # 如果表达量为NaN，认为该酶是未测量的
        if np.isnan(ab):
            non_measured.append(enz)
        else:
            measured[enz] = ab  # 如果表达量有效，记录为已测量的

    # 计算已测量蛋白质的总量
    tot_measured_proteins = np.sum([
        0 if np.isnan(dfomics.loc[ind, condition_id]) else dfomics.loc[ind, condition_id] * dfmws.loc[ind, 'MW']
        for ind in dfomics.index
    ])

    # 计算已测量酶的总量
    tot_measured_enzymes = np.sum([
        measured[ind] * dfmws.loc[ind, 'MW'] for ind in measured.keys()
    ])

    # 计算未测量的蛋白质总量
    tot_non_measured = dftot.loc[condition_id, 'Ptot'] - tot_measured_proteins

    # 计算酶的比例
    prot_pool = tot_non_measured * tot_measured_enzymes / tot_measured_proteins

    # 输出相关结果
    print('tot_proteins         :', dftot.loc[condition_id, 'Ptot'])
    print('tot_measured_proteins:', tot_measured_proteins)
    print('tot_measured_enzymes :', tot_measured_enzymes)
    print('tot_non_measured     :', tot_non_measured)
    print('fraction_enzyme      :', tot_measured_enzymes / tot_measured_proteins)
    print('prot_pool            :', prot_pool)

    # 返回已测量的酶，未测量的酶，以及蛋白质池
    return measured, non_measured, prot_pool


# #基于已测量的酶和蛋白质的丰度比例来估算未测量酶的总丰度。

# In[27]:


def prepare_omics_for_one_condition(dfomics,dftot,dfmws,condition_id,model):
  
    non_measured = []
    measured = {}
    
    assert condition_id in dfomics.columns
    for enz in model.enzymes:
        try: ab = dfomics.loc[enz,condition_id]
        except: ab = np.nan
        
        if np.isnan(ab): non_measured.append(enz)
        else: measured[enz] = ab
    
    tot_measured_proteins = np.sum([0 if np.isnan(dfomics.loc[ind,condition_id]) 
                                    else dfomics.loc[ind,condition_id]*dfmws.loc[ind,'MW'] 
                                    for ind in dfomics.index])
    
    tot_measured_enzymes  = np.sum([measured[ind]*dfmws.loc[ind,'MW'] for ind in measured.keys()])
    tot_non_measured      = dftot.loc[condition_id,'Ptot'] - tot_measured_proteins
    prot_pool = tot_non_measured*tot_measured_enzymes/tot_measured_proteins
    
    print('tot_proteins         :',dftot.loc[condition_id,'Ptot'])
    print('tot_measured_proteins:',tot_measured_proteins)
    print('tot_measured_enzymes :',tot_measured_enzymes)
    print('tot_non_measured     :',tot_non_measured)
    print('fracton_enzyme       :',tot_measured_enzymes/tot_measured_proteins)
    print('prot_pool            :',prot_pool)
    
    return measured, non_measured, prot_pool 


# #基于已测量的酶和蛋白质的丰度比例来估算未测量酶的总丰度。（最后一个参数为 model而不是enzymes，函数从模型中获取酶 ID。）

# In[28]:


measured, non_measured, prot_pool = prepare_omics_for_one_condition(dfomics,dftot,dfmws,'MU',eModel)


# In[29]:


def constrainAbundance(model,measured):

    for prot_id, ab in measured.items():
        rxn_id = 'prot_{0}_exchange'.format(prot_id)
        model.reactions.get_by_id(rxn_id).upper_bound = ab
   
    return model


# #根据给定的蛋白质丰度数据，设置相应蛋白质交换反应的上界为蛋白质的丰度值。

# In[30]:


def constrainPool(model, MWs, measured, non_measured, UB, copy=True):
  
    if copy: 
        model = model.copy()
    
    # Create prot_pool metabolite 
    prot_pool = Metabolite('prot_pool')
    prot_pool.name = prot_pool.id
    
    rxns_to_add  = []
    rxns_to_drop = []
    
    # Loop over enzymes without proteomics data
    for prot in non_measured:
        # Get the corresponding exchange reaction
        prot_exchange_rxn = model.reactions.get_by_id(f'prot_{prot}_exchange')
        
        # Create new draw reaction
        draw_rxn = Reaction(f'draw_prot_{prot}')
        draw_rxn.name = draw_rxn.id
        draw_rxn.add_metabolites({prot_pool: -MWs[prot], list(prot_exchange_rxn.metabolites)[0]: 1})
        
        # Add reactions to respective lists
        rxns_to_add.append(draw_rxn)
        rxns_to_drop.append(prot_exchange_rxn)
    
    # Add new draw reactions to model
    model.add_reactions(rxns_to_add)
    
    # Remove old exchange reactions
    model.remove_reactions(rxns_to_drop)
    
    # Set upper bound of all reactions to np.inf
    for rxn in model.reactions:
        rxn.upper_bound = np.inf
    
    # Create prot_pool_exchange reaction
    rxn_prot_pool_exg = Reaction('prot_pool_exchange')
    rxn_prot_pool_exg.name = rxn_prot_pool_exg.id
    rxn_prot_pool_exg.add_metabolites({prot_pool: 1})
    rxn_prot_pool_exg.lower_bound = 0
    rxn_prot_pool_exg.upper_bound = UB
    
    # Add the prot_pool_exchange reaction to the model (wrap in a list)
    model.add_reactions([rxn_prot_pool_exg])
    
    # Constrain the proteins with measured abundance
    constrainAbundance(model, measured)
    
    return model


# #为未测量的蛋白质设置一个虚拟的蛋白质池，并根据提供的上界来约束这个蛋白质池。

# In[31]:


df_enz_mw = pd.read_csv('data/enzyme_MWs.csv',index_col=0)
df_enz_mw.head(5)


# In[32]:


MWs = {ind:df_enz_mw.loc[ind,'MW'] for ind in df_enz_mw.index}


# In[33]:


ecModel = constrainPool(eModel,MWs, measured, non_measured,prot_pool)


# In[35]:


test_biomass_production(ecModel)

