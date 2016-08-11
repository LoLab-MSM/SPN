''' 
Boolean2Rules

Requires a Boolean model, in the Booleannet format, as input.
Outputs a PySB readable model including header information, Monomers
Initials, Observables, and a set of Boolean rules in mass action format.
'''
from sys import argv
import re
import copy
import pydot
from cStringIO import StringIO
import matplotlib.pyplot as plt
import matplotlib.image as mpimg

class Node(): # define the node class for tree construction
    
    def __init__(self, function, function_nodes, id='', index=1):
        self.id = id
        self.tree = None
        self.index = index
        self.index_name = None
        self.value = None
        self.function = function
        self.function_nodes = function_nodes
        self.true_node = None
        self.false_node = None
        self.true_function = None
        self.false_function = None
        self.reduced_node_set = None
        self.mark = None

    def getTrueNode(self):
        return self.true_node
    def getFalseNode(self):
        return self.false_node
            
    def insertTrue(self, true_function, reduced_node_set, true_id, index):
        self.true_node = Node(true_function, reduced_node_set, true_id, index)
    def insertFalse(self, false_function, reduced_node_set, false_id, index):
        self.false_node = Node(false_function, reduced_node_set, false_id, index)

def constructTree(current_node): # expands the tree via Shannon expansion and computes values of the leaves

    if current_node.index == 1:
        current_node.tree = current_node.id
    nodes = current_node.function_nodes[:]
    id = current_node.id
    index = current_node.index
    function = current_node.function[:]
    function.insert(0, '(')
    function.append(')')
    function_update = []
    true_half = function[:]
    false_half = function[:]
    if nodes != []:
        current_node.index_name = current_node.function_nodes[0]
        for i,token in enumerate(function):
            if token == nodes[0]:
                true_half[i] = 'True'
                false_half[i] = 'False'
        true_half.extend(['and', nodes[0], 'or'])
        false_half.extend(['and', '( not', nodes[0], ')'])
        function_update = true_half + false_half
        true_function = function_update[:]
        false_function = function_update[:]
        for i,bool in enumerate(true_function):
            if bool == nodes[0]:
                true_function[i] = 'True'
        for i,bool in enumerate(false_function):
            if bool == nodes[0]:
                false_function[i] = 'False'                
        current_node.true_function = true_function
        current_node.false_function = false_function         
        current_node.reduced_node_set = nodes[1:]
        current_node.insertTrue(true_function, nodes[1:], id+'1', index+1)
        current_node.insertFalse(false_function, nodes[1:], id+'0', index+1)
        constructTree(current_node.true_node)
        constructTree(current_node.false_node)
    else:
        expression = " ".join(map(str, function))
        current_node.value = eval(expression)
        current_node.index_name = str(eval(expression))

def indexNodes(expansion, indexList=None):  # creates index for use in constructing the ROBDD from the tree
    if indexList == None:
        indexList = {}
    if expansion != None:
        if expansion.index in indexList and expansion not in indexList[expansion.index]:
            indexList[expansion.index].append(expansion)
            indexNodes(expansion.getTrueNode(), indexList)
            indexNodes(expansion.getFalseNode(), indexList)
        else:
            indexList[expansion.index] = []      
            indexList[expansion.index].append(expansion)
            indexNodes(expansion.getTrueNode(), indexList)
            indexNodes(expansion.getFalseNode(), indexList)
    return indexList

def constructBDD(treeRoot): # constructs the ROBDD
    
    BDDroot = copy.deepcopy(treeRoot)
    indexList = indexNodes(BDDroot)
    L = len(indexList)
    
    # define leaf nodes
    True_leaf = Node(None, None, '1', L)
    True_leaf.value = True
    True_leaf.index_name = 'True'  
    False_leaf = Node(None, None, '0', L)
    False_leaf.value = False
    False_leaf.index_name = 'False'

    # set final node level to appropriate leaf node (True/False value)
    for i,j in enumerate(indexList[L-1]):
        if indexList[L-1][i].true_node.value == True:
            indexList[L-1][i].true_node = True_leaf
        else:
            indexList[L-1][i].true_node = False_leaf
        if indexList[L-1][i].false_node.value == True:
            indexList[L-1][i].false_node = True_leaf
        else:
            indexList[L-1][i].false_node = False_leaf  
    n=1
    while L-n != 1:
        
        # find redundant nodes
        for i,j in enumerate(indexList[L-n]):
            if indexList[L-n][i].true_node == indexList[L-n][i].false_node:
                for k,l in enumerate(indexList[L-n-1]):
                    if indexList[L-n-1][k].true_node.id == indexList[L-n][i].id:
                        indexList[L-n-1][k].true_node = indexList[L-n][i].true_node
                        indexList[L-n][i].mark = 'redundant'
                    if indexList[L-n-1][k].false_node.id == indexList[L-n][i].id:
                        indexList[L-n-1][k].false_node = indexList[L-n][i].true_node
                        indexList[L-n][i].mark = 'redundant'
                        
        # remove redundant nodes
        removeList = []
        for i,j in enumerate(indexList[L-n]):
            if indexList[L-n][i].mark == 'redundant':
                removeList.append(indexList[L-n][i])
        for each in removeList:
            indexList[L-n].remove(each)
            
        # find isomorphic nodes
        dup = []
        removeList = []
        for i,j in enumerate(indexList[L-n]):
            if indexList[L-n][i].mark == None:
                indexList[L-n][i].mark = indexList[L-n][i].id
                dup.append(indexList[L-n][i].mark)
                for k,l in enumerate(indexList[L-n]):
                    if (indexList[L-n][k].mark == None) and (indexList[L-n][i].true_node == indexList[L-n][k].true_node) and (indexList[L-n][i].false_node == indexList[L-n][k].false_node):
                        indexList[L-n][k].mark = indexList[L-n][i].id
                        removeList.append(indexList[L-n][k])
        
        # combine isomorphic nodes     
        for each in dup:
            marked = False
            temp_node = None
            for i,j in enumerate(indexList[L-n-1]):
                if indexList[L-n-1][i].true_node.mark == each and marked == False:
                    temp_node = indexList[L-n-1][i].true_node
                    marked = True
                if indexList[L-n-1][i].false_node.mark == each and marked == False:
                    temp_node = indexList[L-n-1][i].false_node
                    marked = True
                                       
                if indexList[L-n-1][i].true_node.mark == each and marked == True:
                    indexList[L-n-1][i].true_node = temp_node
                if indexList[L-n-1][i].false_node.mark == each and marked == True:
                    indexList[L-n-1][i].false_node = temp_node
        for each in removeList:
            indexList[L-n].remove(each)
        n += 1
        
    return BDDroot

def computeTruthTable(function, nodes): # computes a truth table for a Boolean function and a set of nodes
                                        # gives us the initial leaf ordering for FindMinPathOrderHeap
     
    header = copy.deepcopy(nodes)
    header.append('result')
    table = []
    for i in range(2**len(nodes)):
        table.append([])
    for i,node in enumerate(nodes, 1):
        k = len(table)/(2**i)
        count = 1
        value = True
        for j,case in enumerate(table, 1):
            case.append(value)
            if count == k:
                count = 0
                value = not value
            count+=1
         
    for i,each in enumerate(table):
        function_copy = copy.deepcopy(function)
        for j,node in enumerate(nodes):
            for k,token in enumerate(function_copy):
                if token == node:
                    function_copy[k] = each[j]        
        expression = " ".join(map(str, function_copy))
        value = eval(expression)
        table[i].append(value)       
    table.insert(0, header)
    return table

def LeafSwap(nodes, leaves, high, low): # rearranges leaves in accordance to the new node order
    
    l = leaves[:]
    lenN = len(nodes)
    lenL = len(leaves)
    groupSize = 2**(lenN-high+1)
    exchangeSize = 2**(lenN-low)
    ind = 0
    for i in range(lenL/groupSize):
        for j in range(0, groupSize/2, exchangeSize):
            if ind == 0:
                ind = 1
            elif ind == 1:
                ind = 0
                leaves[(i*groupSize + j):(i*groupSize + j)+exchangeSize], leaves[(i*groupSize + j)+groupSize/2-exchangeSize:(i*groupSize + j)+groupSize/2] = \
                leaves[(i*groupSize + j)+groupSize/2-exchangeSize:(i*groupSize + j)+groupSize/2],leaves[(i*groupSize + j):(i*groupSize + j)+exchangeSize]

    return leaves   

def FindMinPathOrderHeap(functions, nodes): # determines a variable order for the minimum number of BDD paths using Heap's algorithm; this is a brute force method 
    
    # initial path reduction count
    newNodes = {}
    for key in functions:
        newNodes[key] = copy.deepcopy(nodes[key])
        N = len(nodes[key])
        leaves = []
        table = computeTruthTable(functions[key], nodes[key])
        for i in range(1, len(table)):
            leaves.append(table[i][len(table[i])-1])
        lenLeaves = len(leaves)
        marks = [1]*(lenLeaves)
        for each in range(N):
            set_size = 2**each
            num_sets = 2**(N-(each))
            counter = 1
            left = None
            right = None       
            for j in range(0, lenLeaves, set_size):
                if counter == 1:
                    left = leaves[j:j+set_size]
                    counter = 2
                else:
                    right = leaves[j:j+set_size]
                    if left == right:
                        for k in range(j, j+set_size):
                            marks[k] = 0
                    counter = 1
        pathcount = 0
        for each in marks:
            pathcount += int(each)
        index = [0 for i in range(N)]
        i = 1
        
        # run through all permutations of node order while counting path reduction
        while i < N:
            if index[i] < i:
                swap = i % 2 * index[i]
                leaves = LeafSwap(nodes[key], leaves, swap+1, i+1)
                marks = [1]*(lenLeaves)
                for each in range(N):
                    set_size = 2**each
                    num_sets = 2**(N-(each))
                    counter = 1
                    left = None
                    right = None       
                    for j in range(0, lenLeaves, set_size):
                        if counter == 1:
                            left = leaves[j:j+set_size]
                            counter = 2
                        else:
                            right = leaves[j:j+set_size]
                            if left == right:
                                for k in range(j, j+set_size):
                                    marks[k] = 0
                            counter = 1
                nodes[key][swap], nodes[key][i] = nodes[key][i], nodes[key][swap]
                paths = 0
                for each in marks:
                    paths += int(each)
                if paths < pathcount:
                    pathcount = paths
                    newNodes[key] = copy.deepcopy(nodes[key])
                    
                index[i] += 1
                i = 1
            else:
                index[i] = 0
                i+= 1
                
    return newNodes

def pathExpansion(expansion, path=None, paths=None):    # lists the paths for a given BDD or tree

    if path == None:
        path = []
    if paths == None:
        paths = []
    if expansion != None:
        if expansion.value != None:
            path.append(expansion.index_name)
            paths.append(path)
        else:
            path.append(expansion.index_name)
            path_t = copy.deepcopy(path)
            path_f = copy.deepcopy(path)
            path_t.append('1')
            pathExpansion(expansion.getTrueNode(), path_t, paths)
            path_f.append('0')
            pathExpansion(expansion.getFalseNode(), path_f, paths)
            
    return paths

def listRules(root):    # lists the rules from the paths in a BDD (or tree)
      
    function_node = root.tree
    paths = pathExpansion(root)
    for i,j in enumerate(paths, 1):
        rule = ''
        print_rule = False
        for k,l in enumerate(j[:-1]):
            if k % 2 == 0:            
                rule += str(l)+'(state = \''
            if k % 2 == 1 and k != len(j)-2:
                rule += str(l)+'\') + '
            if k % 2 == 1 and k == len(j)-2:
                rule += str(l)+'\')'
        if function_node not in j:
            print_rule = True
            if j[len(j)-1] == 'True':
                rule = 'Rule(\''+function_node+str(i)+'\', '+function_node+'(state = \'0\') + '+rule+' >> '+function_node+'(state = \'1\') + '+rule+', on)'
            if j[len(j)-1] == 'False':
                rule = 'Rule(\''+function_node+str(i)+'\', '+function_node+'(state = \'1\') + '+rule+' >> '+function_node+'(state = \'0\') + '+rule+', on)'
        if function_node in j:
            place = j.index(function_node)
            if j[len(j)-1] == 'True' and j[place+1] == '0':
                print_rule = True
                j2 = copy.deepcopy(j)
                j2[place+1] = '1'
                rule2 = ''
                for k,l in enumerate(j2[:-1]):
                    if k % 2 == 0:            
                        rule2 += str(l)+'(state = \''
                    if k % 2 == 1 and k != len(j2)-2:
                        rule2 += str(l)+'\') + '
                    if k % 2 == 1 and k == len(j2)-2:
                        rule2 += str(l)+'\')'
                rule = 'Rule(\''+function_node+str(i)+'\', '+rule+' >> '+rule2+', on)'
            if j[len(j)-1] == 'False' and j[place+1] == '1':
                print_rule = True
                j2 = copy.deepcopy(j)
                j2[place+1] = '0'
                rule2 = ''
                for k,l in enumerate(j2[:-1]):
                    if k % 2 == 0:
                        rule2 += str(l)+'(state = \''
                    if k % 2 == 1 and k != len(j2)-2:
                        rule2 += str(l)+'\') + '
                    if k % 2 == 1 and k == len(j2)-2:
                        rule2 += str(l)+'\')'
                rule = 'Rule(\''+function_node+str(i)+'\', '+rule+' >> '+rule2+', on)'
        if print_rule == True:
            outfile.write(rule+'\n')

def grove(functions, function_nodes): # groups all BDDs (or trees)
    
    tree_list = []
    bdd_list = []
    for keys,values in functions.items():
        tree_root = Node(functions[keys], function_nodes[keys], keys, 1)
        nodeNum = len(function_nodes[keys])
        constructTree(tree_root)
        tree_list.append(tree_root)
        bdd_root = constructBDD(tree_root)
        bdd_list.append(bdd_root)
        
    return bdd_list

def displayIndex(expansion, indexList=None):  # creates index for use in printTree

    if indexList is None:
        indexList = {}
    if expansion != None:
        if expansion.index in indexList:        
            present = False
            for j in indexList[expansion.index]:
                if j == expansion:
                    present = True
            if present == False:
                indexList[expansion.index].append(expansion)
                displayIndex(expansion.getTrueNode(), indexList)
                displayIndex(expansion.getFalseNode(), indexList)
        else:
            indexList[expansion.index] = []                
            present = False
            for j in indexList[expansion.index]:
                if j == expansion:
                    present = True
            if present == False:
                indexList[expansion.index].append(expansion)
                displayIndex(expansion.getTrueNode(), indexList)
                displayIndex(expansion.getFalseNode(), indexList)
                
    return indexList

def printTree(expansion):   # utility for displaying the ROBDD
                            # usage: printTree(constructBDD(Node(function, node_list)))
    
    indexList = displayIndex(expansion)
    nodeList = []
    graph = pydot.Dot(graph_type='digraph')
    for i,j in enumerate(indexList, 1):
        for k,l in enumerate(indexList[i]):
            node = pydot.Node(indexList[i][k].index_name+indexList[i][k].id)
            node.set('label', indexList[i][k].index_name)
            node.set('rank', indexList[i][k].index)
            graph.add_node(node)
        for k,l in enumerate(indexList[i]):
            if i != len(indexList):
                graph.add_edge(pydot.Edge(pydot.Node(indexList[i][k].index_name+indexList[i][k].id), pydot.Node(indexList[i][k].true_node.index_name+indexList[i][k].true_node.id), label='1'))
                graph.add_edge(pydot.Edge(pydot.Node(indexList[i][k].index_name+indexList[i][k].id), pydot.Node(indexList[i][k].false_node.index_name+indexList[i][k].false_node.id), label='0'))

    png = graph.create_png(prog='dot')
    sio = StringIO()
    sio.write(png)
    sio.seek(0)
    img = mpimg.imread(sio)
    imgplot = plt.imshow(img)
    plt.axis('off')
    plt.tight_layout()
    plt.show()

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

orderedNodes = FindMinPathOrderHeap(functions, function_nodes) # minimize the ROBDD paths

# write model header, Monomers, Initials, and Observables
out_name = 'B2R_' + argv[2]
outfile = open(out_name, 'w')
outfile.write('\nfrom pysb import *\n\n')
outfile.write('Model()\n\n')
for key in sorted(initial_states):
    outfile.write('Monomer(\''+key+'\', [\'state\'], {\'state\': [\'0\', \'1\']})\n')
for key in sorted(initial_states):
    init = None
    if initial_states[key] == 'True':
        outfile.write('Initial('+key+'(state = \'0\'), Parameter(\''+key+'0_init\', 0))')
        outfile.write('\n')
        outfile.write('Initial('+key+'(state = \'1\'), Parameter(\''+key+'1_init\', 1))')
        outfile.write('\n')
    if initial_states[key] == 'False':
        outfile.write('Initial('+key+'(state = \'0\'), Parameter(\''+key+'0_init\', 1))')
        outfile.write('\n')
        outfile.write('Initial('+key+'(state = \'1\'), Parameter(\''+key+'1_init\', 0))')
        outfile.write('\n')
outfile.write('\n')
for key in sorted(initial_states):
    outfile.write('Observable(\''+key+'0_obs\', '+key+'(state = \'0\'))')
    outfile.write('\n')
    outfile.write('Observable(\''+key+'1_obs\', '+key+'(state = \'1\'))')
    outfile.write('\n')
outfile.write('\n')
 
# write Rules
BDDs = grove(functions, orderedNodes)   
for each in BDDs:
    listRules(each)
    outfile.write('\n') 
        
outfile.close()  
