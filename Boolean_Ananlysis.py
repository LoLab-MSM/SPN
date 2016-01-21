
from sys import argv
import re
import copy
import datetime
import profile
import random

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

def BoolFuncCalc(node, function, currentState, updateState): ###### Boolean in 1's and 0's
    
    func = function[:]
    for i,each in enumerate(func):
        for item in currentState:
            if each == item[0]:
                func[i] = item[1]
    expression = " ".join(map(str, func))
    value = str(int(eval(expression)))
    for item in updateState:
        if node == item[0]:
            item[1] = value

def FindAttractors():   # synchronous
    
    lenN = len(all_nodes)
    attractors = []
    stateList2 = []
    for i in range(2**lenN):
        stateList2.append([i, None])
    attractor = 0
    k = 0
    while k < len(stateList2):
        if stateList2[k][1] == None:
            stateList2[k][1] = 0
            state = bin(k)[2:].rjust(lenN, '0')
            currentState = [None]*len(initial_states)
            for i,each in enumerate(all_nodes):
                currentState[i] = [all_nodes[i], state[i]]
            updateState = copy.deepcopy(currentState)
            currentRun = []
            currentRun.append(int(state, 2))
            end = False
            while end == False:
                state = ""
                for node,function in functions.items():
                    BoolFuncCalc(node, function, currentState, updateState)
                currentState = copy.deepcopy(updateState)
                for each in currentState:
                    state = state+each[1]
                currentRun.append(int(state, 2))
                if stateList2[int(state, 2)][1] == 0:
                    attractor += 1
                    for every in currentRun:
                        stateList2[every][1] = attractor
                    end = True
                    temp = []
                    for each in currentRun:
                        if each in temp:
                            attractors.append([attractor, currentRun[temp.index(each)+1:]])
                        else:
                            temp.append(each)
                            
                            
                if stateList2[int(state, 2)][1] == None:
                    stateList2[int(state, 2)][1] = 0
                else:
                    for every in currentRun:
                        stateList2[every][1] = stateList2[int(state, 2)][1]
                    end = True
        k += 1   
    print 'Statelist' 
    for each in stateList2:
        print each
    print
    print 'Attractors'
    for each in attractors:
        print each
        
# FindAttractors()

def Synchronous(state):
    
    lenN = len(all_nodes)
    attractor = []
    currentState = [None]*len(initial_states)
    for i,each in enumerate(all_nodes):
        currentState[i] = [all_nodes[i], state[i]]
    updateState = copy.deepcopy(currentState)
    currentRun = []
    currentRun.append(int(state, 2))
    end = False
    while end == False:
        state = ""
        for node,function in functions.items():
            BoolFuncCalc(node, function, currentState, updateState)
        currentState = copy.deepcopy(updateState)
        for each in currentState:
            state = state+each[1]
        currentRun.append(int(state, 2))
        
        temp = []
        for each in currentRun:
            if each in temp:
                attractor.append(currentRun[temp.index(each)+1:])
                end = True
                print currentRun
            else:
                 temp.append(each)

    print 'Attractor'
    for each in attractor:
        print each

# Synchronous('010')

def RandomAsynchronous(state, rounds): # 1 round is each node once

    nodes = all_nodes[:]
    print all_nodes
    
    lenN = len(all_nodes)
    currentState = [None]*len(initial_states)
    for i,each in enumerate(all_nodes):
        currentState[i] = [all_nodes[i], state[i]]
    updateState = copy.deepcopy(currentState)
    
    round = 0
    while round < rounds:
        random.shuffle(nodes)
        print nodes
        
        for node in nodes:
            BoolFuncCalc(node, functions[node], currentState, currentState)
            print currentState

        round += 1

# RandomAsynchronous('010', 10)

def GeneralAsynchronous(state, rounds): # 1 round is one node evaluation
    
    nodes = all_nodes[:]
    print all_nodes
    
    lenN = len(all_nodes)
    currentState = [None]*len(initial_states)
    for i,each in enumerate(all_nodes):
        currentState[i] = [all_nodes[i], state[i]]
    updateState = copy.deepcopy(currentState)
    
    round = 0
    while round < rounds:
        node = random.choice(nodes)
        print node
         
        BoolFuncCalc(node, functions[node], currentState, currentState)
        print currentState

        round += 1

GeneralAsynchronous('011', 10)







