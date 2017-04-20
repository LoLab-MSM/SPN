
from sys import argv
import re
import copy
import datetime
import profile
import random
import operator

# read and tokenize the Boolean model
Boolean_model = open(argv[2], 'r')

all_nodes = []          # list of all nodes in system
initial_states = {}     # initial states for each node
function_nodes = {}     # list of nodes incident on each node (Boolean function variables)         
functions = {}          # tokenized list of all nodes, operations, and parentheses in each function

in_rules = False
for line in Boolean_model: # checks for 'rules' text block between triple quotes
    if re.match('.*"""', line,) and in_rules == False:
        in_rules = True
    if re.match('"""', line,) and in_rules == True:
        in_rules = False

    if in_rules == True: # extract tokens for each Boolean function
        if re.match('.*\*.*', line):
            line = line.replace('(', '( ')
            line = line.replace(')', ' )')
            for dNode in all_nodes:
                if re.match('.*' + str(dNode) + '\*.*', line):
                    node_list = []
                    token_list = []
                    fuction_list = []
                    for iNode in all_nodes:
                        if re.match('.*\s*' + str(iNode) + '(\s.*|$)', line):
                            node_list.append(iNode)
                        token_list = re.split(' ', line.rstrip('\n'))
                    function_nodes[dNode] = node_list
                    index_of_eq = 0
                    for index, token in enumerate(token_list):
                        if token == '=':
                            index_of_eq = index
                    function_list = token_list[index_of_eq+1:]
                    functions[dNode] = function_list
        else:       # extract initial states
            if line.rstrip('\n')[-4:] == 'True' or line.rstrip('\n')[-5:] == 'False' or line.rstrip('\n')[-6:] == 'Random':
                line_list = re.split('\s*=\s*', line.rstrip('\n'))
                initial = line_list[-1]
                line_list = line_list[:-1]
                for node in line_list:
                    initial_states[node] = initial
                    all_nodes.append(node)

def BoolFuncCalc(node, function, currentState): ###### Boolean in 1's and 0's
    
    func = function[:]
    for i,each in enumerate(func):
        for item in currentState:
            if each == item[0]:
                func[i] = item[1]
    expression = " ".join(map(str, func))
    value = str(int(eval(expression)))
    for item in currentState:
        if node == item[0]:
            item[1] = value

def GeneralAsynchronous(state, rounds): # 1 round is one node evaluation
    
    nodes = all_nodes[:]
    trace = []
    
    lenN = len(all_nodes)
    currentState = [None]*len(initial_states)
    for i,each in enumerate(all_nodes):
        currentState[i] = [all_nodes[i], str(int(eval(initial_states[each])))]
    trace.append(copy.deepcopy(currentState))
    round = 0
    while round < rounds:
        node = random.choice(nodes)
        BoolFuncCalc(node, functions[node], currentState)
        trace.append(copy.deepcopy(currentState))

        round += 1
    
    return trace

output = GeneralAsynchronous(initial_states, 10)

for each in output:
    print each

# def FindAttractors():   # synchronous; no memory reduction
#     
#     lenN = len(all_nodes)
#     attractors = []
#     statelist2 = []
#     for i in range(2**lenN):
#         statelist2.append(None)
#     attractor = 0
#     k = 0
#     while k < len(statelist2):
#         if statelist2[k] == None:
#             statelist2[k] = 0
#             state = bin(k)[2:].rjust(lenN, '0')
#             currentState = [None]*len(initial_states)
#             updateState = [None]*len(initial_states)
#             for i,each in enumerate(all_nodes):
#                 currentState[i] = [all_nodes[i], state[i]]
#                 updateState[i] = [all_nodes[i], state[i]]
# #             updateState = copy.deepcopy(currentState)
#             currentRun = []
#             currentRun.append(int(state, 2))
#             end = False
#             while end == False:
#                 state = ""
#                 for node,function in functions.items():
#                     BoolFuncCalc(node, function, currentState, updateState)
# #                 currentState = copy.deepcopy(updateState)
#                 currentState[0][:],currentState[1][:],currentState[2][:] = updateState[0][:],updateState[1][:],updateState[2][:]
#                 for each in currentState:
#                     state = state+each[1]
#                 currentRun.append(int(state, 2))
#                 if statelist2[int(state, 2)] == 0:
#                     attractor += 1
#                     for every in currentRun:
#                         statelist2[every] = attractor
#                     end = True
#                     temp = []
#                     for each in currentRun:
#                         if each in temp:
#                             attractors.append([attractor, 0, currentRun[temp.index(each)+1:]])
#                         else:
#                             temp.append(each)
#                             
#                             
#                 if statelist2[int(state, 2)] == None:
#                     statelist2[int(state, 2)] = 0
#                 else:
#                     for every in currentRun:
#                         statelist2[every] = statelist2[int(state, 2)]
#                     end = True
#         attractors[statelist2[k]-1][1] += 1
#         k += 1
#         
#     print 'Statelist' 
#     for i,each in enumerate(statelist2):
#         print i, each
#     print
#     
#     print 'Attractors'
#     for each in attractors:
#         print each
#
# # FindAttractors()
# 
# print 
# print '-----------------------------' 
# print
# 
# def FindAttractorsCompressed():   # synchronous; with memory reduction
#     
#     lenN = len(all_nodes)
#     stateN = 2**lenN
#     statelist = []
#     compressedlist = []
#     templist = {}
#     attractors = []
#     attractor = 0    
#     k = 0
#     s = 0
#     while k < stateN:
#         statelist.append([k, None, templist.get(k)])
#         if statelist[k+s][2] == None:
#             statelist[k+s][2] = 0
#             state = bin(k)[2:].rjust(lenN, '0')
#             currentState = [None]*len(initial_states)
#             updateState = [None]*len(initial_states)
#             for i,each in enumerate(all_nodes):
#                 currentState[i] = [all_nodes[i], state[i]]
#                 updateState[i] = [all_nodes[i], state[i]]
# #             updateState = copy.deepcopy(currentState)
#             currentRun = []
#             currentRun.append(int(state, 2))
#             end = False
#             while end == False:
#                 state = ""
#                 for node,function in functions.items():
#                     BoolFuncCalc(node, function, currentState, updateState)
# #                 currentState = copy.deepcopy(updateState)
#                 currentState[0][:],currentState[1][:],currentState[2][:] = updateState[0][:],updateState[1][:],updateState[2][:]
#                 for each in currentState:
#                     state = state+each[1]
#                 currentRun.append(int(state, 2))
#                 if int(state, 2) <= k:
#                     if statelist[int(state, 2)+s][2] == 0:
#                         attractor += 1
#                         for every in currentRun:
#                             if every <= k:
#                                 statelist[every+s][2] = attractor
#                             else:
#                                 templist[every] = attractor
#                         end = True
#                         temp = []
#                         for each in currentRun:
#                             if each in temp:
#                                 attractors.append([attractor, 0, currentRun[temp.index(each)+1:]])
#                             else:
#                                 temp.append(each)
#                     if statelist[int(state, 2)+s][2] == None:
#                         statelist[int(state, 2)][2] = 0
#                     else:
#                         for every in currentRun:
#                             if every <= k:
#                                 statelist[every+s][2] = statelist[int(state, 2)+s][2]
#                             else:
#                                 templist[every] = statelist[int(state, 2)][2]
#                         end = True
#                 else:
#                     if templist.get(int(state, 2)) == 0:
#                         attractor += 1
#                         for every in currentRun:
#                             if every <= k:
#                                 statelist[every+s][2] = attractor
#                             else:
#                                 templist[every] = attractor
#                         end = True
#                         temp = []
#                         for each in currentRun:
#                             if each in temp:
#                                 attractors.append([attractor, 0, currentRun[temp.index(each)+1:]])
#                             else:
#                                 temp.append(each)
#                     if templist.get(int(state, 2)) == None:
#                         templist[int(state, 2)] = 0
#                     else:
#                         for every in currentRun:
#                             if every <= k:
#                                 statelist[every+s][2] = templist[int(state, 2)]
#                             else:
#                                 templist[every] = templist[int(state, 2)]
#                         end = True
#                         
#             if k > 0 and statelist[k+s][2] == statelist[k-1+s][2]:
#                 statelist[k-1+s][1] = statelist[k+s][0]
#                 del statelist[k+s]
#                 s -= 1
#         else:
#             del templist[k]
#             if k > 0 and statelist[k+s][2] == statelist[k-1+s][2]:
#                 statelist[k-1+s][1] = statelist[k+s][0]
#                 del statelist[k+s]
#                 s -= 1
#                 
#         attractors[statelist[k+s][2]-1][1] += 1    
#         k += 1
#     for i in range(len(statelist)):
#         if statelist[i][1] == None:
#             statelist[i][1] = statelist[i][0]
#             
#     print 'Statelist' 
#     for each in statelist:
#         print each
# 
#     print
# #     for each in templist:
# #         print each, templist[each]
#     
#     print 'Attractors'
#     for each in attractors:
#         print each
# 
# # FindAttractorsCompressed()
# 
# def Synchronous(state):
#     
#     lenN = len(all_nodes)
#     attractor = []
#     currentState = [None]*len(initial_states)
#     for i,each in enumerate(all_nodes):
#         currentState[i] = [all_nodes[i], state[i]]
#     updateState = copy.deepcopy(currentState)
#     currentRun = []
#     currentRun.append(int(state, 2))
#     end = False
#     while end == False:
#         state = ""
#         for node,function in functions.items():
#             BoolFuncCalc(node, function, currentState)
#         currentState = copy.deepcopy(updateState)
#         for each in currentState:
#             state = state+each[1]
#         currentRun.append(int(state, 2))
#         
#         temp = []
#         for each in currentRun:
#             if each in temp:
#                 attractor.append(currentRun[temp.index(each)+1:])
#                 end = True
#                 print currentRun
#             else:
#                  temp.append(each)
# 
#     print 'Attractor'
#     for each in attractor:
#         print each
# 
# # Synchronous('010')
# 
# def RandomAsynchronous(state, rounds): # 1 round is each node once
# 
#     nodes = all_nodes[:]
#     print all_nodes
#     
#     lenN = len(all_nodes)
#     currentState = [None]*len(initial_states)
#     for i,each in enumerate(all_nodes):
#         currentState[i] = [all_nodes[i], state[i]]
#     updateState = copy.deepcopy(currentState)
#     
#     round = 0
#     while round < rounds:
#         random.shuffle(nodes)
#         print nodes
#         for node in nodes:
#             BoolFuncCalc(node, functions[node], currentState)
#             print currentState
# 
#         round += 1
# 
# # RandomAsynchronous('001', 10)




