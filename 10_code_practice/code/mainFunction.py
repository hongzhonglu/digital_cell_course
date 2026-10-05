# -*- coding: utf-8 -*-
# -*- python 3 -*-
# -*- hongzhong Lu -*-

# note: this file is not latest, will be not updated further.
# 2021.12.02




# Import packages
import numpy as np
import pandas as pd


'''general function for easy use of python'''
def splitAndCombine(gene, rxn, sep0, moveDuplicate=False):
    ## one rxn has several genes, this function was used to splite the genes
    ## used for the dataframe data

    gene = gene.fillna('NA')  # fill the NaN with 'NA'
    gene0 = gene.tolist()
    rxn0 = rxn.tolist()
    s1 = list()
    s2 = list()
    for i in range(len(gene0)):
        s1 = s1 + [rxn0[i]] * len(gene0[i].split(sep0))
        s2 = s2 + gene0[i].split(sep0)
    df0 = pd.DataFrame({'V1': s1,
                        'V2': s2}
                       )
    if moveDuplicate == True:
        df00 = df0.drop_duplicates()
    else:
        df00 = df0
    return df00


def getSimilarTarget(rxn_yeast0,rxn_newGPR0,ss):
    from fuzzywuzzy import fuzz
    from fuzzywuzzy import process
    rxn_yeast1 = np.array(rxn_yeast0)  # np.ndarray()
    rxn_yeast2 = rxn_yeast1.tolist()
    rxn_yeast3 = pd.Series((v[0] for v in rxn_yeast2))
    rxn_newGPR1 = np.array(rxn_newGPR0)  # np.ndarray()
    rxn_newGPR2 = rxn_newGPR1.tolist()
    rxn_newGPR3 = pd.Series((v[0] for v in rxn_newGPR2))
    similarTarget = [None] * ss
    for i in range(ss):
        similarTarget[i] = process.extract(rxn_newGPR3[i], rxn_yeast3, limit=2)

    return similarTarget
'''
#example
newMet = pd.read_excel('new metabolite for check.xlsx')
newMet0 = newMet[['name_unify']]
gemMet = pd.read_excel('unique metabolite in yeastGEM.xlsx')
gemMet0 = gemMet[['Description_simple']]
ss0 = len(newMet0)
similarTarget0 = getSimilarTarget(gemMet0,newMet0,ss=ss0)
'''


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
'''
w=['a','b','c']
v=[1,2,3]
s=[3,1,2,4]
singleMapping(w,v,s,dataframe=False)
'''


def multiMapping (description, item1, item2, dataframe=True, sep=";", removeDuplicates=True):
    """get multiple description of from item1 for item2 based on mapping"""
    #description = w
    #item1 = v
    #item2 = testData
    #used for the list data
    if dataframe:
        description = description.tolist()
        item1 = item1.tolist()
        item2 = item2.tolist()
    else:
        pass
    result = [None]*len(item2)
    for i in range(len(item2)):
        if item2[i] in item1:
            index0 = [description[index] for index in range(len(item1)) if item1[index] == item2[i]]
            if removeDuplicates:
                index1 = pd.unique(index0).tolist()
            else:
                index1 = index0
            result[i] = sep.join(str(e) for e in index1) #string cat
        else:
            result[i] = None
    return result

'''
# example data to test all the above function
df1 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'C' : ['foo', 'foo', 'foo', 'bar', 'bar', 'bar'] * 2}
                   )

df2 = pd.DataFrame({'A' : ['one', 'one', 'two', 'three'] * 3,
                    'B' : ['A', 'B', 'C'] * 4,
                    'D' : np.random.randn(12)})


df2['C'] = singleMapping(df1['C'], df1['A'], df2['A'])
df2['C'] = multiMapping(df1['C'], df1['A'], df2['A'])
'''



def produceGeneList(model0):
  #produce the gene list from GEM
  genelist = []
  for i in model0.genes:
      print(i)
      genelist.append(i.id)
  return genelist


def saveExcel(infile, outfile):
    '''
    function to save the dataframe into xlsx format
    :param infile:
    :param outfile:
    :return:
    '''
    writer = pd.ExcelWriter(outfile)
    infile.to_excel(writer,'Sheet1')
    writer.save()

