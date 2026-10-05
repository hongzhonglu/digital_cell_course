#example to reconstruct model based on cobra
from __future__ import print_function
from cobra import Model, Reaction, Metabolite
from cobra.sampling import sample


"""Building a Model"""
#creat model and reaction
model = Model('example_model')

#define the metabolite
def metabolite_input(short_name, formula0, detail_name, compartment):
    return Metabolite(
    short_name,
    formula= formula0,
    name=detail_name,
    compartment=compartment)

GLU = metabolite_input('GLU',None,None,'c')
GLY = metabolite_input('GLY',None,None,'c')
ACE = metabolite_input('ACE',None,None,'c')
ETH = metabolite_input('ETH',None,None,'c')
CO2= metabolite_input('CO2',None,None,'c')
NADH = metabolite_input('NADH',None,None,'c')
ATP = metabolite_input('ATP',None,None,'c')
BIO = metabolite_input('BIO',None,None,'c')


#define the reaction
def reaction_input(reactionID, composition):
    reaction = Reaction(reactionID)
    reaction.name = reactionID
    #Adding metabolites to a reaction
    reaction.add_metabolites(composition)
    #add the reaction into model
    model.add_reactions([reaction])
    return model

r1 = 'glucose_input'
dict1 = {GLU: 1}

r2 = 'acetaldehyde_sythesis'
dict2 = {GLU: -1.5,
         ACE:1,
         CO2:0.5,
         NADH:0.5,
         ATP:0.5
}

r3 = 'ethanol_synthesis'
dict3 = {GLU:-1,
         NADH:-0.5,
         ETH: 1}

r4 = 'glycerol_synthesis'
dict4 = {GLU:-1,
         NADH: -0.333,
         ATP:-0.333,
         GLY:1}

r5 = 'biomass_synthesis'
dict5 = {GLU:-1.12,
         ATP: -2.42,
         BIO:1,
         CO2:0.12,
         NADH:0.15}

r6 = 'biomass_export'
dict6 = {BIO:-1}

r7 = 'glycerol_export'
dict7 = {GLY:-1}

r8 = 'ethanol_export'
dict8 = {ETH:-1}

r9 = 'acetaldehyde_export'
dict9 = {ACE:-1}

r10 = 'co2_export'
dict10 = {CO2:-1}

rxnID = [r1,r2,r3,r4,r5,r6,r7,r8,r9,r10]
rxnFormula0 = [dict1,dict2,dict3, dict4, dict5,dict6,dict7,dict8,dict9,dict10]


for i,x in enumerate(rxnID):
    print(i, x)
    model = reaction_input(x, rxnFormula0[i])

'''add related info to reactions
    gpr
    subsystem
    lower bound
    upper bound'''
gpr = ['','a','b','c', '(m and n)','','','','','']
subsystem = ['Exchange','','','','','','','','','']
upper_bound = [1000]*10

for i in range(len(gpr)):
    print(i)
    model.reactions[i].gene_reaction_rule = gpr[i]
    model.reactions[i].subsystem = subsystem[i]
    model.reactions[i].upper_bound = upper_bound[i]

'''medium'''
medium = model.medium
medium["glucose_input"] = 10
model.medium = medium


'''objective function'''
model.solver = 'glpk'
#model.solver = 'cplex'
model.objective = 'biomass_synthesis'
solution = model.optimize()
solution.status
solution.fluxes
print(model.summary())




# Now there are things in the model
print('%i reaction' % len(model.reactions))
print('%i metabolites' % len(model.metabolites))
print('%i genes' % len(model.genes))


#We can iterate through the model objects to observe the contents
# Iterate through the the objects in the model
print("Reactions")
print("---------")
for x in model.reactions:
    #print("%s : %s" % (x.id, x.reaction))
    if x.id == "glycerol_synthesis":
        print(x.reaction)


print("")
print("Metabolites")
print("-----------")
for x in model.metabolites:
    print('%9s : %s' % (x.id, x.formula))

print("")
print("Genes")
print("-----")
for x in model.genes:
    associated_ids = (i.id for i in x.reactions)
    print("%s is associated with reactions: %s" %
          (x.id, "{" + ", ".join(associated_ids) + "}"))


'''Advanced function'''
#sampling
s = sample(model, 100)
s = sample(model, 100, method="achr")
s.head()
s.columns

# visualization
import seaborn as sns
import matplotlib.pyplot as plt
fig = plt.figure(figsize=(4,3))
sns.distplot(s['ethanol_export'])
sns.distplot(s['glucose_input'])
fig.legend(labels=['ethanol_export','glucose_input'])
plt.show()

from cobra.sampling import OptGPSampler, ACHRSampler
achr = ACHRSampler(model, thinning=10)
optgp = OptGPSampler(model, processes=4)
s1 = achr.sample(100)
s2 = optgp.sample(100)

fig = plt.figure(figsize=(4,3))
sns.distplot(s1['ethanol_export'])
sns.distplot(s2['ethanol_export'])
fig.legend(labels=['achr','optgp'])
plt.show()




co = model.problem.Constraint(model.reactions.glucose_input.flux_expression, lb=5)
model.add_cons_vars([co])
optgp = OptGPSampler(model, processes=4)
s3 = optgp.sample(100)

fig = plt.figure(figsize=(4,3))
sns.distplot(s2['ethanol_export'])
sns.distplot(s3['ethanol_export'])
fig.legend(labels=['no constraint','with constraint'])
plt.show()

fig = plt.figure(figsize=(4,3))
sns.distplot(s2['glucose_input'])
sns.distplot(s3['glucose_input'])
fig.legend(labels=['no constraint','with constraint'])
plt.show()


