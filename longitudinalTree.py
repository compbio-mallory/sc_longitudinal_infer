import pandas as pd
import json
import argparse
from collections import OrderedDict
import itertools
from itertools import combinations, chain

class Tree:
    def __init__(self):
        self.node = None
        self.data = None

''' Convert input matrix to pandas. '''
def convertToPandas(tsv_fileName):
    pd_fromFile = pd.read_csv(tsv_fileName, sep='\t', index_col=0)
    return pd_fromFile

''' Obtain subsets using a python module itertools. '''
def find_subsets(set_1, snvs, noOfsnv):
    finalList = []
    #subsetList = list(map(set, itertools.combinations(snvs, noOfsnv)))
    subsetList = itertools.combinations(snvs, noOfsnv)
    for subset1 in subsetList:
        #print("subset1 ",subset1)
        if (set_1.issubset(subset1)): # Only add list of subsets from previous time.
            finalList.append(subset1)
    return finalList

''' This method validates the obtained subclone by checking subset and superset of the obtained subclones. '''
def validate_subclones(temp_sc, dict_For_tree):
    count_usc = 0
    temp_dict = dict_For_tree.copy()
    for key in temp_sc:
        key_value = set(temp_sc[key])
        #sets_not_to_include = 
        #print(key_value)
        for input_key in dict_For_tree:
            if key in input_key:
                #print("Here is the key ..",key)
                value = set(dict_For_tree[input_key])
                #print(value)
                for sets in key_value:
                    if (set(sets) != value and set(sets).issuperset(value)): # Check if the obtained subclone is a superset of the subclone in previous step.
                        #if(sets_present != [] and "False" in sets_present
                        count_usc = count_usc + 1
                        #print(key," ",sets," --> ",value)
                        temp_dict[key+"_"+"usc"+str(count_usc)] = list(sets)
                    else:
                        #temp_temp_dict = temp_dict.copy()
                        for temp_key in temp_dict: # This else condition helps in checking for all subclones in a given time to make sure the unobserved subclones remains proper superset and subset.
                            #print(temp_key)
                            if key in temp_key and "usc" in temp_key:
                                #print(value," -- ",temp_dict[temp_key])
                                if(value == set(temp_dict[temp_key])):
                                    #print(value," -- ",temp_dict[temp_key])
                                    temp_temp_dict = dict(temp_dict)
                                    del temp_temp_dict[temp_key]
                                    temp_dict = temp_temp_dict.copy() # Update the dict with the proper superclones and subclones. Remove any subclone failing this condition.
    return temp_dict

''' Find the unobserved subclones. '''
def find_unobserved_subclones(dict_For_tree):
    temp_unobserved_sc = {} # will save all subsets here and separate them later.
    for key_1 in dict_For_tree: # iterating through input_tree
        time_point_k = key_1.split('_')[1] # we just want to get the times and not the whole key (ex. t_1_c1)
        set_1 =	set(dict_For_tree[key_1])
        #subset = set()
        #unobserved_sc = set()
        #print(key_1," length of subset ",len(subset))
        #print(" =========================== ")
        for key_2 in dict_For_tree: # comparing subclones of two times pairwise to find subsets.
            if key_1 == key_2:
                continue
            time_point_k1 = key_2.split('_')[1]
            if int(time_point_k)+1 == int(time_point_k1): # condition to check only between t_1 and t_2 and not between t_1 and t_3. Maintaining the time stepwise.
                #print(len(subset))
                set_2 = set(dict_For_tree[key_2])
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
            #if row[0] == 1:
            #    print(row[0])
    return dict_For_tree

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (Input to the Longitudinal tree)")
args = parser.parse_args()
tsv_input = args.input

inputToTree = convertToPandas(tsv_input)
dict_For_tree = convert_InputToDict(inputToTree)
print("Step 1 to get the clones arranged ... ",dict_For_tree)
unobserved_sc = find_unobserved_subclones(dict_For_tree)

#print(new_dict_tree)
#set_1 = {'pos1'}
#set_2 = {'pos2', 'pos1', 'pos4'}
#if(set_1 != set_2 and set_1.issubset(set_2)):
#    print("Proper subset")
