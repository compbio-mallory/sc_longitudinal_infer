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
        #print(key_value)
        for input_key in dict_For_tree: # Comparing subclones with the original dictionary to validate the subset and superset of clones.
            if key in input_key:
                #print("Here is the key ..",key)
                value = set(dict_For_tree[input_key])
                #print(value)
                for sets in key_value:
                    if (set(sets) != value and set(sets).issuperset(value)): # Check if the obtained subclone is a superset of the subclone in previous step.
                        count_usc = count_usc + 1
                        #print(key," ",sets," --> ",value)
                        temp_dict[key+"_"+"usc"+str(count_usc)] = list(sets)
                    else:
                        for temp_key in temp_dict: # This else condition helps in checking for all subclones in a given time to make sure the unobserved subclones remains proper superset and subset.
                            #print(temp_key)
                            if key in temp_key and "usc" in temp_key:
                                #print(value," -- ",temp_dict[temp_key])
                                temp_temp_dict = dict(temp_dict)
                                del temp_temp_dict[temp_key]
                                temp_dict = temp_temp_dict.copy() # Update the dict with the proper superclones and subclones. Remove any subclone failing this condition.
    return temp_dict

''' Find the unobserved subclones. '''
def find_unobserved_subclones(dict_For_tree):
    temp_unobserved_sc = {} # will save all subsets here and separate them later with the function validate_subclones()
    for key1 in dict_For_tree: # iterating through input_tree. Remember keys here are combination of time and cluster.
        time_point_k = key1.split('_')[1] # we just want to get the times and not the whole key (ex. t_1_c1)
        set_1 =	set(dict_For_tree[key1]) # set_1 indicates the SNVs in a given time and cluster 
        #print(key_1," length of subset ",len(subset))
        #print(" =========================== ")
        for key2 in dict_For_tree: # comparing subclones of two times pairwise to find subsets.
            if key1 == key2:
                continue
            time_point_k1 = key2.split('_')[1]
            if int(time_point_k)+1 == int(time_point_k1): # condition to check only between t_1 and t_2 and not between t_1 and t_3. Maintaining the time stepwise.
                #print(len(subset))
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
            

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (Input to the Longitudinal tree)")
args = parser.parse_args()
tsv_input = args.input

inputToTree = convertToPandas(tsv_input)
clones_in_time = convert_InputToDict(inputToTree)
print("Step 1 to get the clones arranged ... ",clones_in_time)
unobserved_sc_tree = find_unobserved_subclones(clones_in_time)
tree_with_children = construct_tree(unobserved_sc_tree)
print(" Tree with children ",tree_with_children)
edge_weight_dict = calculate_edge_weight(unobserved_sc_tree, tree_with_children)
print(" Tree with edge weights ",edge_weight_dict)

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
