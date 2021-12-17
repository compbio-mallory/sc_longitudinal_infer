import pandas as pd
import json
import argparse
from collections import OrderedDict
import itertools
from itertools import combinations, chain

#class Tree:
#    def __init__(self,node,data, tree_dict,children='None'):
#        self.tree_dict = tree_dict
#        self.node = node
#        self.data = data
#        self.children = [] # Add the nodes as children for reference
#        if children is not None:
#            for child in children:
#                self.add_child(child)
        #if(children != 'None' or children != []):
        #    self.tree_dict[self.node] = self.children
        #return self.tree_dict
        #print("Node ",self.node," Children ",self.children)
#        self.create_tree(node, children, tree_dict)

#    def __repr__(self):
#        return self.tree_dict
    
#    def add_child(self, node):
        #assert isinstance(node, Tree)
#        self.children.append(node)

#    def create_tree(self, node, children, tree_dict):
        #def __setitem__(self, key, value):
        #setattr(self, key, value)
#        if(children is not None or children != []):
            #setattr(self.tree_dict, self.node, self.children)
#            tree_dict[node] : children

''' Convert input matrix to pandas. '''
def convertToPandas(tsv_fileName):
    pd_fromFile = pd.read_csv(tsv_fileName, sep='\t', index_col=0)
    return pd_fromFile

''' Obtain subsets using a python module itertools. '''
def find_subsets(t1_set, t2_set, noOfsnv): # set_1 are the SNVs in t. snvs are the ones in t+1. noOfsnv is the subset size.
    finalList = []
    #subsetList = list(map(set, itertools.combinations(snvs, noOfsnv)))
    subsetList = itertools.combinations(t2_set, noOfsnv)
    for subset1 in subsetList:
        #print("subset1 ",subset1)
        if (t1_set != subset1 and t1_set.issubset(subset1)): # Only add list of proper subsets from previous time.
            finalList.append(subset1)
    return finalList

''' This method validates the obtained subclone by checking subset and superset of the obtained subclones. '''
def validate_subclones(temp_sc, dict_For_tree):
    count_usc = 0
    temp_dict = dict_For_tree.copy()
    for key in temp_sc: # Remember keys here are time and clusters. E.g. t1_c1
        key_value = set(temp_sc[key]) # key_value indicates the SNVs present.
        #print("SNVs for a key ",key_value)
        for input_key in dict_For_tree: # Comparing subclones with the original dictionary to validate the subset and superset of clones.
            if key in input_key:
                #print("Here is the key ..",key)
                value = set(dict_For_tree[input_key])
                #print(" SNVs of original tree ",value)
                for sets in key_value:
                    if (set(sets) != value and set(sets).issuperset(value)): # Check if the obtained subclone is a superset of the subclone in previous step.
                        count_usc = count_usc + 1
                        #print("Checking if subclone is superset of time t ",key," ",sets," --> ",value)
                        temp_dict[key+"_"+"usc"+str(count_usc)] = list(sets)
                        #print(temp_dict)
                        #print("-------------------")
                    #else:
                    #    for temp_key in temp_dict: # This else condition helps in checking for all subclones in a given time to make sure the unobserved subclones remains proper superset and subset.
                            #print(temp_key)
                    #        if key in temp_key and "usc" in temp_key:
                    #            print(" Rechecking the subset superset ",value," -- ",temp_dict[temp_key])
                    #            temp_temp_dict = dict(temp_dict)
                    #            del temp_temp_dict[temp_key]
                    #            temp_dict = temp_temp_dict.copy() # Update the dict with the proper superclones and subclones. Remove any subclone failing this condition.
                    #            print("**********",temp_dict)
    return temp_dict

''' Find the unobserved subclones. '''
def find_unobserved_subclones(dict_For_tree):
    temp_unobserved_sc = {} # will save all subsets here and separate them later with the function validate_subclones()
    for key1 in dict_For_tree: # iterating through input_tree. Remember keys here are combination of time and cluster.
        time_point_k = key1.split('_')[1] # we just want to get the times and not the whole key (ex. t_1_c1)
        set_1 =	set(dict_For_tree[key1]) # set_1 indicates the SNVs in a given time and cluster 
        #print("key 1 ",key1," length of subset ",len(set_1))
        #print(" =========================== ")
        for key2 in dict_For_tree: # comparing subclones of two times pairwise to find subsets.
            if key1 == key2:
                continue
            time_point_k1 = key2.split('_')[1]
            if int(time_point_k)+1 == int(time_point_k1): # condition to check only between t_1 and t_2 and not between t_1 and t_3. Maintaining the time stepwise.
                #print("Key 2 ",key2)
                set_2 = set(dict_For_tree[key2]) # set_2 indicates the SNVs in t+1 and each cluster in t+1 
                if(set_1 != set_2 and set_1.issubset(set_2)): #checking proper subsets
                    usc = find_subsets(set_1,set_2, len(set_2)-1) # find all the subsets for a matching subclone
                    temp_unobserved_sc_key = 't_'+time_point_k
                    if temp_unobserved_sc_key not in temp_unobserved_sc: # This if else statement is to add list of subsets in a particular time. (ex. subsets between t_1 and t_2) 
                        temp_unobserved_sc[temp_unobserved_sc_key] = usc
                    else:
                        temp_sc = temp_unobserved_sc[temp_unobserved_sc_key]
                        temp_sc.extend(usc)
                        temp_unobserved_sc[temp_unobserved_sc_key] = temp_sc
            else:
                continue
            
    #print(temp_unobserved_sc)
    unobserved_sc = validate_subclones(temp_unobserved_sc,dict_For_tree) # check if the subsets obtained are proper subsets and supersets.
    print("Dict with the unobserved clones ",unobserved_sc)
    return unobserved_sc

''' Method to create a dictionary from the input matrix which will help in creating tree. '''
def convert_InputToDict(input_df):
    #for cluster, row in GprimeDf.iterrows():
    dict_For_tree = OrderedDict()
    input_df_dict = input_df.to_dict('records')
    colNames = input_df_dict[0].keys()
    clusterNo = 0
    for row in input_df.itertuples():
        #print(row.Index)
        clusterNo = clusterNo+1
        for col in colNames:
            #print(col)
            dictList = []
            if getattr(row,col) == 1:
                dictKey = row.Index+"_c"+str(clusterNo)
                if dictKey in dict_For_tree:
                    tempList = dict_For_tree[dictKey]
                    tempList.append(col)
                    dict_For_tree[dictKey] = tempList
                else:
                    dictList.append(col)
                    dict_For_tree[dictKey] = dictList
            #print(dict_For_tree)
            #print(input_df.loc[row.Index, col])
    return dict_For_tree

''' Construct a tree as dictionary where nodes are the keys and their corresponding children nodes are values. '''
def construct_tree(unobserved_sc_tree):
    #Compare keys of each pairs and check if their values are subsets.
    #If yes then add them as children. Adding children is basically adding edges.
    tree_dict = {}
    for k_time in unobserved_sc_tree:
        k_time_value = set(unobserved_sc_tree[k_time])
        time_point_k = k_time.split('_')[1]
        k_time_children = []
        #k_time_node = Tree(k_time,unobserved_sc_tree[k_time])
        for k1_time in unobserved_sc_tree:
            if k_time == k1_time:
                continue
            time_point_k1 = k1_time.split('_')[1]
            #k1_time_node = Tree(k1_time, unobserved_sc_tree[k1_time],tree_dict)
            if int(time_point_k)+1 == int(time_point_k1): #checking neighboring time points
                k1_time_value = set(unobserved_sc_tree[k1_time])
                if(k_time_value.issubset(k1_time_value)):
                    k_time_children.append(k1_time)
            elif(int(time_point_k) == int(time_point_k1) and 'usc' in k1_time): #checking unobserved subclones here
                k1_time_value = set(unobserved_sc_tree[k1_time])
                if(k_time_value.issubset(k1_time_value)):
                    k_time_children.append(k1_time)
            else:
                continue
        tree_dict[k_time] = k_time_children
        #tree_dict = Tree(k_time,unobserved_sc_tree[k_time],tree_dict,k_time_children)
    return tree_dict

''' For each child node check if there are new SNVs and update the edge weight if new SNVs are found. ''' 
def calculate_edge_weight(unobserved_sc_tree, tree_with_children):
    edge_weight_dict = {}
    for parent_node in tree_with_children:
        children = tree_with_children[parent_node]
        for child in children:
            parent_snvCount = len(unobserved_sc_tree[parent_node])
            child_snvCount = len(unobserved_sc_tree[child])
            keytuple = (parent_node, child)
            if(parent_snvCount == child_snvCount):
                edge_weight_dict[keytuple] = 0
            elif(child_snvCount > 0):
                edge_weight_dict[keytuple] = child_snvCount - parent_snvCount
    return edge_weight_dict

''' Selecting the node s with maximum number of connections.'''
def select_nodeS(node_list,tree_with_children,nodes_to_remove):
    maxValue = 0
    s_node = ""
    if len(node_list) != 0:
        for node in node_list:
            if nodes_to_remove != [] and node in nodes_to_remove:
                continue
            node_t1_conn = len(tree_with_children[node])
            if maxValue == 0 or maxValue < node_t1_conn:
                maxValue = node_t1_conn
                s_node = node
    return s_node

''' Listing the nodes to remove from the tree. '''
def del_nodes(s_node, node_list):
    nodes_to_remove = []
    for node in node_list:
        if s_node != node:
            nodes_to_remove.append(node)
    return nodes_to_remove

''' Selecting A_bar  '''
def select_A_bar(tree_with_children, t1, s_node):
    A_bar_nodes = []
    s_node_children = tree_with_children[s_node]
    for node in tree_with_children:
        if t1 in node:
            if node not in s_node_children and 'usc' not in node:
                A_bar_nodes.append(node)
        else:
            continue
    return A_bar_nodes

''' Finding any unobserved subclones between time t and union of s_node and A_bar or only s_node. '''
def enroll_unobservedSc(s_node, A_bar, tree_with_children, unobserved_sc_tree):
    new_t1_nodes = [] # union of s_node and A_bar. Or just s_node when A_bar is empty.
    new_t1_nodes_dict = {}
    t_nodes_dict = {}
    nodes_to_add = {}

    print(" S node is ",s_node," A bar nodes are ",A_bar)
    if len(A_bar) != 0:
        new_t1_nodes.append(s_node)
        new_t1_nodes.extend(A_bar)
    else:
         new_t1_nodes.append(s_node)
         
    for node in new_t1_nodes:
        new_t1_nodes_dict[node] = unobserved_sc_tree[node]

    t_k = s_node.split('_')[1]
    for node in unobserved_sc_tree: #Taking only the nodes in time t
        if "t_"+t_k in node and 'usc' not in node:
            t_nodes_dict[node] = unobserved_sc_tree[node]

    usc_count = 0
    for t_node in t_nodes_dict:
        t_node_snv = set(t_nodes_dict[t_node])
        for t1_node in new_t1_nodes_dict:
            #print("S node and A_bar nodes .. ",t1_node)
            t1_node_snv = set(new_t1_nodes_dict[t1_node])
            if t_node_snv != t1_node_snv and t_node_snv.issubset(t1_node_snv):
                usc = find_subsets(t_node_snv,t1_node_snv, len(t1_node_snv)-1)
                for sc in usc:
                    usc_count = usc_count+1
                    key = s_node+"_usc"+str(usc_count)
                    #print("Nodes to add key ",key) #checking below that the node is only proper subset and proper superset
                    if t_node_snv != set(sc) and t_node_snv.issubset(set(sc)) and t1_node_snv != set(sc) and t1_node_snv.issuperset(set(sc)):
                        nodes_to_add[key] = sc
    
    #print("New t1 nodes ",new_t1_nodes)
    #print(" t1 dict ",new_t1_nodes_dict)
    #print(" t dict ",t_nodes_dict)
    #print(nodes_to_add)
    return nodes_to_add
        
''' Greedy algorithm to remove/add unobserved clones and preserve mutations in the forward direction. '''
def forward_mutation_one(tree_with_children,unobserved_sc_tree):
    nodes_to_remove = []
    nodes_to_add_dict = {}
    for node_t in tree_with_children:
        t_children = tree_with_children[node_t]
        #print(t_children)
        usc_list = []
        for subclones in t_children:
            if 'usc' in subclones:
                usc_list.append(subclones)
        #print(node_t+" contains "+str(usc_list))
        s_node = select_nodeS(usc_list, tree_with_children, nodes_to_remove)
        #print(s_node) # Check if this is empty before further operations.
        if del_nodes(s_node, usc_list) != []:
            nodes_to_remove.extend(del_nodes(s_node, usc_list))
        print("Nodes to remove ",set(nodes_to_remove)) # Check always if a node is in this list before further operations.
        
        t_k = node_t.split('_')[1]
        for node_t1 in tree_with_children:
            if node_t == node_t1:
                continue
            t1_k1 = node_t1.split('_')[1]
            if int(t_k)+1 == int(t1_k1):
                if (s_node != ""):
                    A_bar_nodes = select_A_bar(tree_with_children, "t_"+t1_k1, s_node)
                    print("A bar nodes ",A_bar_nodes)
                    nodes_to_add = enroll_unobservedSc(s_node, A_bar_nodes, tree_with_children, unobserved_sc_tree)
                    nodes_to_add_dict.update(nodes_to_add)
                    #print("Nodes to add ",nodes_to_add)
            t1_children = tree_with_children[node_t1]
    return nodes_to_remove, nodes_to_add_dict

''' This part searches for unconnected nodes in t+1 (u) and connect them to a node in t (v).'''
def forward_mutation_two(tree_with_children, unobserved_sc_tree):
    # The following lines searches if any node in t+1 is unconnected to any node in t.
    check_if_node_exists_list = []
    children_list = [] # This keeps count if node in t+1 is present as children of any nodes in t. 
    for node_t in unobserved_sc_tree:
        t_children = tree_with_children[node_t]
        if t_children == []:
            continue
        #print(" Node t child ",t_children)
        t_k = node_t.split('_')[1]
        for node_t1 in unobserved_sc_tree:
            if node_t == node_t1:
                continue
            t1_k1 = node_t1.split('_')[1]
            if int(t_k)+1 == int(t1_k1):
                #print(" Node t1 ",node_t1)
                if len(check_if_node_exists_list) != 0:
                    for node in check_if_node_exists_list:
                        if node in t_children:
                            check_if_node_exists_list.remove(node)
                if node_t1 not in t_children and node_t1 not in children_list:
                    check_if_node_exists_list.append(node_t1)
                if node_t1 in t_children:
                    children_list.append(node_t1)
    #print(" If child exists node ",set(check_if_node_exists_list))

    if check_if_node_exists_list == []: # Nothing else to do if all nodes are already connected.
        return None

    node_u_v = {}
    for lone_node in check_if_node_exists_list: # Adding node u and node v in this loop prior to checking subsets.
        t1_k1 = lone_node.split('_')[1]
        max_mutations = 0
        t1_snv_set = set(unobserved_sc_tree[lone_node])
        for t_node in unobserved_sc_tree:
            if t_node == lone_node:
                continue
            t_node_children = tree_with_children[t_node]
            if 'usc' in t_node or lone_node in t_node_children:
                continue
            t_k = t_node.split('_')[1]
            if int(t_k)+1 == int(t1_k1):
                t_snv_set = set(unobserved_sc_tree[t_node])
                if t1_snv_set != t_snv_set and t_snv_set.issubset(t1_snv_set):
                    snv_count = len(t_snv_set)
                    if max_mutations == 0:
                        max_mutations = snv_count
                        node_u_v[lone_node] = t_node
                    elif snv_count >= max_mutations:
                        max_mutations = snv_count
                        node_u_v[lone_node] = t_node

    print(" Node_u_v ",node_u_v)
    return node_u_v

''' This method evaluates the cluster to check if any mutation is missed at a time point and then later added at another time.
This is a clustering error which we want to detect.'''
def evaluate_clustering(tree_with_children, unobserved_sc_tree):
    cluster_mistake = 0
    print(" Evaluating cluster =============== ")
    for node_t in unobserved_sc_tree:
        if 'usc' in node_t: # only consider the observed subclones in this time 
            continue
        t_node_snv = unobserved_sc_tree[node_t] # note the mutations in time t
        t_node_children = tree_with_children[node_t] # checking if the nodes are connected
        t_k = node_t.split('_')[1]
        s = 2 # this variable will help to check in time point >=2 of current time point
        for node_t1 in unobserved_sc_tree:
            if node_t == node_t1:
                continue
            t1_k1 = node_t1.split('_')[1]
            if int(t_k)+s > int(t1_k1):
                continue
            #print(t1_k1+" "+t_k+" "+str(s))
            if int(t1_k1) == int(t_k)+s: #checking in s>=2 time point
                #print(" Node t1 ",node_t1)
                t1_node_snv = unobserved_sc_tree[node_t1]
                #print(" t node SNV ",t_node_snv," t1 node SNV ",t1_node_snv)
                if (t_node_children == [] and set(t_node_snv) == set(t1_node_snv)): # the t node should be unconnected and similar to another node in another time point
                    cluster_mistake = cluster_mistake+1
                    #print("Cluster mistake ",cluster_mistake)
        s = s+1
    return cluster_mistake
    #print(" Evaluating cluster ",cluster_mistake)
    
''' This method looks for back mutations, connects the nodes(x,y) and records the edges. '''
def back_mutation(tree_with_children, unobserved_sc_tree):
    # The following lines searches if any node in t+1 is unconnected to any node in t.
    print("Back mutation steps =======================")
    check_if_node_exists_list = []
    children_list = []
    for node_t in unobserved_sc_tree:
        t_children = tree_with_children[node_t]
        if t_children == []:
            continue
        #print(" Node t child ",t_children)
        t_k = node_t.split('_')[1]
        for node_t1 in unobserved_sc_tree:
            if node_t == node_t1:
                continue
            t1_k1 = node_t1.split('_')[1]
            if int(t_k)+1 == int(t1_k1):
                #print(" Node t1 ",node_t1) 
                if len(check_if_node_exists_list) != 0:
                    for node in check_if_node_exists_list:
                        if node in t_children:
                            check_if_node_exists_list.remove(node)
                if node_t1 not in t_children and node_t1 not in children_list:
                    check_if_node_exists_list.append(node_t1)
                if node_t1 in t_children:
                    children_list.append(node_t1)
    print(check_if_node_exists_list)
    if check_if_node_exists_list == []: # Nothing else to do if all nodes are already connected. 
        return None

    node_x_y = {}
    check_if_node_exists_set = set(check_if_node_exists_list)
    for lone_node in check_if_node_exists_set: # Adding node x and node y in this loop prior to checking supersets.
        t1_k1 = lone_node.split('_')[1]
        min_mutations = 0
        print(" Lone node in t+1 ", lone_node)
        t1_snv_set = set(unobserved_sc_tree[lone_node])
        for t_node in unobserved_sc_tree:
            if t_node == lone_node:
                continue
            t_node_children = tree_with_children[t_node]
            if 'usc' in t_node or lone_node in t_node_children:
                continue
            t_k = t_node.split('_')[1]
            if int(t_k)+1 == int(t1_k1):
                print(" t node : ",t_node," its children ",t_node_children)
                t_snv_set = set(unobserved_sc_tree[t_node])
                if t1_snv_set != t_snv_set and t_snv_set.issuperset(t1_snv_set):
                    snv_count = len(t_snv_set) # need to check the snv count here
                    print(" SNV count ",snv_count)
                    if min_mutations == 0:
                        min_mutations = snv_count # Initially first value are minimum
                        print(" Min mutations ",min_mutations)
                        node_x_y[lone_node] = t_node # keeping the same key so that new connection gets updated.
                    elif snv_count < min_mutations:
                        min_mutations = snv_count # Choosing the minimum value
                        print(" Min mutations ",min_mutations)
                        node_x_y[lone_node] = t_node
                    elif snv_count == min_mutations:
                        min_mutations = snv_count # if two same values exist then choose randomly
                        print(" Min mutations ",min_mutations)
                        node_x_y[lone_node] = t_node

    print(" Node_x_y ",node_x_y)
    return node_x_y
    
''' Method to update the longitudinal tree with the nodes. ''' 
def update_longitudinal_tree(unobserved_sc_tree, nodes_to_remove, nodes_to_add_dict, u_v_conn, x_y_conn):
    # Update unobserved_sc_tree by removing nodes_to_remove and adding values from nodes_to_add_dict.
    # Then call construct_tree function to get the tree with updated nodes.
    # Then add u_v_conn nodes to the constructed tree.
    temp_sc_tree = unobserved_sc_tree.copy()
    print("Nodes to del ",nodes_to_remove)
    print(" Nodes to add ",nodes_to_add_dict)
    for node in nodes_to_remove:
        if node in temp_sc_tree:
            print("Node to del ",node)
            del temp_sc_tree[node]

    for key in nodes_to_add_dict:
        temp_sc_tree[key] = nodes_to_add_dict[key]

    print(" After adding and deleting node ",temp_sc_tree)
    updated_tree = construct_tree(temp_sc_tree)
    if (u_v_conn == {} or u_v_conn == None) and (x_y_conn == {} or x_y_conn == None):
        return updated_tree
    else:
        # check for u_v_conn and x_y_conn
        updated_tree_copy = updated_tree.copy()
        if(u_v_conn != {} or u_v_conn != None):
            for u_v in u_v_conn:
                u_v_key = u_v_conn[u_v]
                updated_tree_copy[u_v_key] = [u_v]
        if(x_y_conn != {} or x_y_conn != None):
            for x_y in x_y_conn:
                print("Updating Back mutation tree ",x_y)
                t_node_key = x_y_conn[x_y]
                updated_tree_copy[t_node_key] = [x_y]
        return updated_tree_copy
    
parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (Input to the Longitudinal tree)")
parser.add_argument("-inputType", "--inputType",dest="inputType", help="Input type can be matrix or dict")
args = parser.parse_args()

input_type = args.inputType
if input_type == 'matrix':
    tsv_input = args.input
    inputToTree = convertToPandas(tsv_input) # accept input as both matrix and dictionary. Prefer using dictionary.
    clones_in_time = convert_InputToDict(inputToTree)
elif input_type == 'dict':
    json_input = args.input
    with open(json_input) as json_file:
        dict_input = json.load(json_file)
    clones_in_time = OrderedDict(dict_input)
    
print("Step 1 to get the clones arranged ... ",clones_in_time)
unobserved_sc_tree = find_unobserved_subclones(clones_in_time)

tree_with_children = construct_tree(unobserved_sc_tree)
print(" Tree with children ",tree_with_children)
#edge_weight_dict = calculate_edge_weight(unobserved_sc_tree, tree_with_children)
#print(" Tree with edge weights ",edge_weight_dict)
nodes_to_remove, nodes_to_add_dict = forward_mutation_one(tree_with_children, unobserved_sc_tree)
u_v_conn = forward_mutation_two(tree_with_children, unobserved_sc_tree)

# the method for counting clustering algo mistake before back mutation goes here.
cluster_mistake = evaluate_clustering(tree_with_children, unobserved_sc_tree)

x_y_conn = back_mutation(tree_with_children, unobserved_sc_tree) # consider this as the back mutation edges collection to analyze later
updated_tree = update_longitudinal_tree(unobserved_sc_tree, nodes_to_remove, nodes_to_add_dict, u_v_conn, x_y_conn)
print("Updated Longitudinal tree ",updated_tree)

#left = Tree("left","left_data")
#middle = Tree("middle","middle_data")
#right = Tree("right","right_data")
#root = Tree("root","root_data", [left, middle, right])
#root.children = [left, middle, right]
#print(root)

#print(new_dict_tree)
#set_1 = {'pos1'}
#set_2 = {'pos2', 'pos1', 'pos4'}
#if(set_1 != set_2 and set_1.issubset(set_2)):
#    print("Proper subset")
