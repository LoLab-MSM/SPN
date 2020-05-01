
# Boolean SPN model from Albert and Othmer 2002

from boolean2 import Model
import numpy as np
import matplotlib.pyplot as plt

rules = """
SLP_1 = False
wg_mRNA_1 = False
WG_1 = False
en_mRNA_1 = True
EN_1 = False
hh_mRNA_1 = True
HH_1 = False
ptc_mRNA_1 = False
PTC_1 = False
PH_1 = False
SMO_1 = False
ci_mRNA_1 = False
CI_1 = False
CIA_1 = False
CIR_1 = False
 
SLP_2 = False
wg_mRNA_2 = False
WG_2 = False
en_mRNA_2 = False
EN_2 = False
hh_mRNA_2 = False
HH_2 = False
ptc_mRNA_2 = True
PTC_2 = False
PH_2 = False
SMO_2 = False
ci_mRNA_2 = True
CI_2 = False
CIA_2 = False
CIR_2 = False
 
SLP_3 = True
wg_mRNA_3 = False
WG_3 = False
en_mRNA_3 = False
EN_3 = False
hh_mRNA_3 = False
HH_3 = False
ptc_mRNA_3 = True
PTC_3 = False
PH_3 = False
SMO_3 = False
ci_mRNA_3 = True
CI_3 = False
CIA_3 = False
CIR_3 = False
 
SLP_4 = True
wg_mRNA_4 = True
WG_4 = False
en_mRNA_4 = False
EN_4 = False
hh_mRNA_4 = False
HH_4 = False
ptc_mRNA_4 = True
PTC_4 = False
PH_4 = False
SMO_4 = False
ci_mRNA_4 = True
CI_4 = False
CIA_4 = False
CIR_4 = False

SLP_1* = SLP_1
wg_mRNA_1* = (CIA_1 and SLP_1 and not CIR_1) or (wg_mRNA_1 and (CIA_1 or SLP_1) and not CIR_1)
WG_1* = wg_mRNA_1
en_mRNA_1* = (WG_4 or WG_2) and not SLP_1
EN_1* = en_mRNA_1
hh_mRNA_1* = EN_1 and not CIR_1
HH_1* = hh_mRNA_1
ptc_mRNA_1* = CIA_1 and not EN_1 and not CIR_1
PTC_1* = ptc_mRNA_1 or (PTC_1 and not HH_4 and not HH_2)
PH_1* = (ptc_mRNA_1 or (PTC_1 and not HH_4 and not HH_2)) and (hh_mRNA_4 or hh_mRNA_2)
SMO_1* = not (ptc_mRNA_1 or (PTC_1 and not HH_4 and not HH_2)) or hh_mRNA_4 or hh_mRNA_2
ci_mRNA_1* = not EN_1
CI_1* = ci_mRNA_1
CIA_1* = CI_1 and (SMO_1 or hh_mRNA_4 or hh_mRNA_2)
CIR_1* = CI_1 and not SMO_1 and not hh_mRNA_4 and not hh_mRNA_2

SLP_2* = SLP_2
wg_mRNA_2* = (CIA_2 and SLP_2 and not CIR_2) or (wg_mRNA_2 and (CIA_2 or SLP_2) and not CIR_2)
WG_2* = wg_mRNA_2
en_mRNA_2* = (WG_1 or WG_3) and not SLP_2
EN_2* = en_mRNA_2
hh_mRNA_2* = EN_2 and not CIR_2
HH_2* = hh_mRNA_2
ptc_mRNA_2* = CIA_2 and not EN_2 and not CIR_2
PTC_2* = ptc_mRNA_2 or (PTC_2 and not HH_1 and not HH_3)
PH_2* = (ptc_mRNA_2 or (PTC_2 and not HH_1 and not HH_3)) and (hh_mRNA_1 or hh_mRNA_3)
SMO_2* = not (ptc_mRNA_2 or (PTC_2 and not HH_1 and not HH_3)) or hh_mRNA_1 or hh_mRNA_3
ci_mRNA_2* = not EN_2
CI_2* = ci_mRNA_2
CIA_2* = CI_2 and (SMO_2 or hh_mRNA_1 or hh_mRNA_3)
CIR_2* = CI_2 and not SMO_2 and not hh_mRNA_1 and not hh_mRNA_3

SLP_3* = SLP_3
wg_mRNA_3* = (CIA_3 and SLP_3 and not CIR_3) or (wg_mRNA_3 and (CIA_3 or SLP_3) and not CIR_3)
WG_3* = wg_mRNA_3
en_mRNA_3* = (WG_2 or WG_4) and not SLP_3
EN_3* = en_mRNA_3
hh_mRNA_3* = EN_3 and not CIR_3
HH_3* = hh_mRNA_3
ptc_mRNA_3* = CIA_3 and not EN_3 and not CIR_3
PTC_3* = ptc_mRNA_3 or (PTC_3 and not HH_2 and not HH_4)
PH_3* = (ptc_mRNA_3 or (PTC_3 and not HH_2 and not HH_4)) and (hh_mRNA_2 or hh_mRNA_4)
SMO_3* = not (ptc_mRNA_3 or (PTC_3 and not HH_2 and not HH_4)) or hh_mRNA_2 or hh_mRNA_4
ci_mRNA_3* = not EN_3
CI_3* = ci_mRNA_3
CIA_3* = CI_3 and (SMO_3 or hh_mRNA_2 or hh_mRNA_4)
CIR_3* = CI_3 and not SMO_3 and not hh_mRNA_2 and not hh_mRNA_4

SLP_4* = SLP_4
wg_mRNA_4* = (CIA_4 and SLP_4 and not CIR_4) or (wg_mRNA_4 and (CIA_4 or SLP_4) and not CIR_4)
WG_4* = wg_mRNA_4
en_mRNA_4* = (WG_3 or WG_1) and not SLP_4
EN_4* = en_mRNA_4
hh_mRNA_4* = EN_4 and not CIR_4
HH_4* = hh_mRNA_4
ptc_mRNA_4* = CIA_4 and not EN_4 and not CIR_4
PTC_4* = ptc_mRNA_4 or (PTC_4 and not HH_3 and not HH_1)
PH_4* = (ptc_mRNA_4 or (PTC_4 and not HH_3 and not HH_1)) and (hh_mRNA_3 or hh_mRNA_1)
SMO_4* = not (ptc_mRNA_4 or (PTC_4 and not HH_3 and not HH_1)) or hh_mRNA_3 or hh_mRNA_1
ci_mRNA_4* = not EN_4
CI_4* = ci_mRNA_4
CIA_4* = CI_4 and (SMO_4 or hh_mRNA_3 or hh_mRNA_1)
CIR_4* = CI_4 and not SMO_4 and not hh_mRNA_3 and not hh_mRNA_1
"""

model = Model(rules, mode='sync')
model.initialize()
model.iterate(steps=30)

#------------------------------------------------------- for node in model.data:
    #---------------------------------------------- print node, model.data[node]

#cycle = model.report_cycles()
nodeList = []
cell_1_states = []
cell_2_states = []
cell_3_states = []
cell_4_states = []
stateList = str(model.states[model.detect_cycles()[0]])
stateList = stateList.split(' ', 1)[1]
stateList = stateList.split(', ')
for i,state in enumerate(stateList):
    if i % 4 == 0:
        nodeList.append(state.rsplit('_', 1)[0])
        cell_1_states.append(state.rsplit('=', 1)[1])
    if i % 4 == 1:
        cell_2_states.append(state.rsplit('=', 1)[1])
    if i % 4 == 2:
        cell_3_states.append(state.rsplit('=', 1)[1])       
    if i % 4 == 3:
        cell_4_states.append(state.rsplit('=', 1)[1])        


column_labels = list('1234')
row_labels = nodeList

string = ''

for node in cell_1_states:
    if node == 'True':
        string = string + '1'
    elif node =='False':
        string = string + '0'
    string = string + ' '
string = string[:-1]
string = string + '; '
for node in cell_2_states:
    if node == 'True':
        string = string + '1'
    elif node =='False':
        string = string + '0'
    string = string + ' '
string = string[:-1]
string = string + '; '
for node in cell_3_states:
    if node == 'True':
        string = string + '1'
    elif node =='False':
        string = string + '0'
    string = string + ' '
string = string[:-1]
string = string + '; '
for node in cell_4_states:
    if node == 'True':
        string = string + '1'
    elif node =='False':
        string = string + '0'
    string = string + ' '
string = string[:-1]

data = np.matrix(string).T
data = np.array(data)

row_labels = np.array(row_labels)

############
# print np.where(row_labels == 'wg_mRNA')[0][0]
# ordered_row_labels = ['wg_mRNA', 'WG', 'en_mRNA', 'EN', 'hh_mRNA', 'HH', 'ptc_mRNA', 'PTC', 'PH', 'SMO', 'ci_mRNA', 'CI', 'CIA', 'CIR']
ordered_row_labels = ['en_mRNA', 'EN', 'wg_mRNA', 'WG', 'ptc_mRNA', 'PTC', 'ci_mRNA', 'CI', 'CIR', 'hh_mRNA', 'HH', 'PH', 'SMO', 'CIA']
# ordered_row_labels = ['CIA', 'CIR', 'CI', 'EN', 'HH', 'PH', 'PTC', 'SMO', 'WG', 'ci_mRNA', 'en_mRNA', 'hh_mRNA', 'ptc_mRNA', 'wg_mRNA']
ordered_data = []
for label in ordered_row_labels:
    index = np.where(row_labels == label)[0][0]
    ordered_data.append(data[index])
ordered_data = np.array(ordered_data)
############

fig, ax = plt.subplots()
# heatmap = ax.pcolor(data, cmap=plt.cm.Blues, edgecolors='black', linewidths=.2)
heatmap = ax.pcolor(ordered_data, cmap=plt.cm.Blues, edgecolors='black', linewidths=.2)

# ax.set_xticks(np.arange(data.shape[1]) + 0.5, minor = False)
# ax.set_yticks(np.arange(data.shape[0]) + 0.5, minor = False)
ax.set_xticks(np.arange(ordered_data.shape[1]) + 0.5, minor = False)
ax.set_yticks(np.arange(ordered_data.shape[0]) + 0.5, minor = False)

# plt.ylim(0,len(row_labels)
plt.ylim(0,len(ordered_row_labels))
ax.invert_yaxis()
ax.xaxis.tick_top()

ax.set_xticklabels(column_labels, minor=False)
# ax.set_yticklabels(row_labels, minor=False)
# ax.set_yticklabels(ordered_row_labels, minor=False)
ax.set_yticklabels([label.split('_')[0] for label in ordered_row_labels], minor=False)

plt.show()
