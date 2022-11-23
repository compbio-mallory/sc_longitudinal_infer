"""
This file is to siimulate the input for longitudinal tree.
"""
import pandas as pd
import json
import argparse
from collections import defaultdict

def json_to_dict(json_file):
    with open(json_file) as j_file:
        jsonDict = json.load(j_file)
    return jsonDict

def convertToPandas(fileName):
    pd_fromFile = pd.read_csv(fileName, compression='gzip', sep='\t', index_col=0)
    return pd_fromFile

def chromPosToIntDict(chromPosList,initial_matrix):
    chromPosToInt_dict = {}
    #print(len(chromPosList))
    i = 1
    for pos in chromPosList:
        chromPosToInt_dict[pos] = i
        i=i+1
    return chromPosToInt_dict
    
# Add a dictionary method here for modifying the pos values. This is just to compare the values and test. Replace pos with a int number. 
''' Function to build the initial matrix from the output of SCG. This will look similar to matrix D. '''
def final_initial_matrix(gp_table, cellToClusterDict):
    initial_matrix = pd.DataFrame()
    chromPosList = []
    for cell in cellToClusterDict: # key is the Cell ID so we check for each cell ID its SNV types on a specific position                                                                                   
        clusterID = cellToClusterDict[cell]
        for row in gp_table.itertuples():
            if clusterID == str(row.Index) and getattr(row,'event_type') == 'snv' and getattr(row,'probability') > 0.80:
                pos = getattr(row,'event_id')
                chromPosList.append(pos)
                SNV = getattr(row,'event_value')
                initial_matrix.loc[cell,pos] = SNV
        #for index, row in gp_table.iterrows():                                                                                                                                                            
        #    if clusterID == str(index) and row['event_type'] == 'snv': # Choose the cluster the cell ID belongs to and get the position(event_id) and SNV(event_value) to fill the matrix                 
        #        pos = row['event_id']                                                                                                                                                                     
        #        SNV = row['event_value']                                                                                                                                                                  
                #print("Event value .. ",event_value)                                                                                                                                                      
        #        initial_matrix.loc[cell,pos] = SNV
    print(" Set values for chrom pos ",set(chromPosList))
    chromPosToInt_dict = chromPosToIntDict(set(chromPosList),initial_matrix)
    initial_matrix.rename(columns=chromPosToInt_dict, inplace=True) # rename chrom pos with simple numbers for testing.
    return initial_matrix.fillna(0),chromPosToInt_dict

''' returns the voted SNV with the values 0,1,2'''
def voting_algo(cell_SNV_list):
    count1 = 0
    count0 = 0
    for snv in cell_SNV_list:
        if (snv == 1 or snv == 2):
            count1 = count1 + 1
        elif(snv == 0):
            count0 = count0 + 1

    if (count0 == count1 or count1 > count0):
        return 1
    else:
        return 0
    
''' Function to get the G' matrix '''
def build_GprimeMatrix(gp_table, final_matrix, chromPosToInt_dict, cellToClusterDict):
    clusterToCellDict = {} # Building a dictionary with cluster IDs as keys and cell IDs as correspondig values                                                                                            
    for cellID in cellToClusterDict:
        clusterID = cellToClusterDict[cellID]
        if clusterID in clusterToCellDict:
            temp_list = clusterToCellDict[clusterID];
            temp_list.append(cellID)
        else:
            clusterToCellDict[clusterID] = [cellID]
    #print(clusterToCellDict.keys())                                                                                                                                                                       
    #posSet = set(gp_table['event_id']) # This was before converting the chrom pos to numbers.
    posSet = final_matrix.columns.values.tolist()
    GprimeDf = pd.DataFrame()
    for cluster in clusterToCellDict:
        for row in gp_table.itertuples():
            if cluster == str(row.Index) and getattr(row,'event_type') == 'snv' and getattr(row,'probability') > 0.80:
                #prob = getattr(row,'probability')
                #if math.ceil(prob) == 1:
                    pos = getattr(row,'event_id')
                    SNV = getattr(row,'event_value')
                    if SNV == 2.0:
                        SNV = 1.0
                    GprimeDf.loc[cluster,pos] = SNV
    GprimeDf.rename(columns=chromPosToInt_dict, inplace=True) # rename chrom pos with simple numbers for testing.
    return GprimeDf, posSet

''' Getting all the SNVs of cells.'''
def cell_pos_SNV(initial_matrix, posSet):
    cell_pos_dict = {}
    for row in initial_matrix.itertuples():
        cellId = row.Index
        pos_list = []
        for pos in posSet:
            #if ':' not in pos:
            #    continue
            #print(initial_matrix.loc[row.Index, pos])
            if initial_matrix.loc[row.Index, pos] == 1.0 or initial_matrix.loc[row.Index, pos] == 2.0:
                pos_list.append(pos)
                cell_pos_dict[cellId] = pos_list
    return cell_pos_dict

''' Getting only the positions with SNVs. '''
def filter_cell_pos_SNV(cell_pos_dict, GprimeDf, clusterToCellDict, posSet):
    filtered_cell_pos_dict = {}
    for clusterId in clusterToCellDict:
        cells_list = clusterToCellDict[clusterId]
        for cell in cell_pos_dict:
            if cell in cells_list:
                pos_list = []
                for pos in cell_pos_dict[cell]:
                    if GprimeDf.loc[clusterId, pos] == 1.0:
                        pos_list.append(pos)
                        filtered_cell_pos_dict[cell] = pos_list
    return filtered_cell_pos_dict

# The way count_sc is included is wrong. Need to fix this. Test this whole algorithm by our simple example. Trace back from the input tree in Figure 1 to create a matrix D.
def simulated_input_tree(clusterToCellDict, timeToCellIdDict, filtered_cell_pos_dict):
    time_sc_cell = {}
    count_sc = 0
    for clusterId in clusterToCellDict:
        cell_list = clusterToCellDict[clusterId]
        count_sc = count_sc+1
        for cell in cell_list:
            #count_sc = count_sc+1
            for time in timeToCellIdDict:
                #count_sc = count_sc+1
                key="t_"+str(time)+"_sc"+str(count_sc)
                time_cell_list = timeToCellIdDict[time]
                subclone_list = []
                if cell in time_cell_list:
                    subclone_list.append(cell)
                    time_sc_cell[key] = subclone_list
                    
    #print(time_sc)
    time_sc_snv = {}
    for time_sc in time_sc_cell:
        time_sc_cell_list = time_sc_cell[time_sc]
        time_snv_list = []
        for cellId in time_sc_cell_list:
            if cellId in filtered_cell_pos_dict.keys():
                cellId_snv = filtered_cell_pos_dict[cellId]
                if time_sc in time_sc_snv:
                    temp_snv_list = time_sc_snv[time_sc]
                    temp_snv_list.extend(cellId_snv)
                    time_sc_snv[time_sc] = temp_snv_list
                else:
                    time_snv_list.extend(cellId_snv)
                    time_sc_snv[time_sc] = time_snv_list
    print(time_sc_snv)
    return time_sc_snv

    #simulated_tree_matrix = pd.DataFrame() # If at any time we want the input to be a matrix then use this. Otherwise dictionary is faster and better to use.
    #for time in time_sc_snv:
    #    snv_list = time_sc_snv[time]
    #    for snv in snv_list:
    #        simulated_tree_matrix.loc[time,snv] = 1
    #print(simulated_tree_matrix.fillna(0))
    #return simulated_tree_matrix.fillna(0)
    
parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (similar to input to SCG)")
parser.add_argument("-gp", "--genotype",dest ="genotype", help="Genotype posterior from SCG")
parser.add_argument("-cellCluster", "--cellCluster",dest ="cellCluster", help="Cell cluster from SCG")
parser.add_argument("-timepoints", "--timepoints",dest ="timepoints", help="Cell time points")
parser.add_argument("-output", "--output",dest ="output", help="Output file to save")
args = parser.parse_args()
#gp_path = input("Enter the path to genotype_posteriors.tsv from SCG: ") # Path now is ../no_doublet/genotype_posteriors.tsv                                                                               
gp_path = args.genotype
gp_table = convertToPandas(gp_path) # gp_table is the dataframe with values from genotype_posteriors.tsv

print(gp_table)
#cellToClusterDict = json_to_dict('cellToCluster.json')
#timeToCellIdDict = json_to_dict('cell_timepoints.json')
cellToClusterDict = json_to_dict(args.cellCluster)
timeToCellIdDict = json_to_dict(args.timepoints)

clusterToCellDict = {}
for cellId in cellToClusterDict:
    clusterID = cellToClusterDict[cellId]
    if clusterID in clusterToCellDict:
        cell_list = clusterToCellDict[clusterID]
        cell_list.append(cellId)
        clusterToCellDict[clusterID] = cell_list
    else:
        temp_cell_list = []
        temp_cell_list.append(cellId)
        clusterToCellDict[clusterID] = temp_cell_list

#print(clusterToCellDict)

initial_matrix, chromPosToInt_dict = final_initial_matrix(gp_table, cellToClusterDict) # Populated matrix similar to G''                                                                                                              
print("=================== INITIAL MATRIX ========================")
print(initial_matrix)
initial_matrix.to_csv('initial_matrix.tsv',sep='\t')

GprimeDf, posSet = build_GprimeMatrix(gp_table, initial_matrix, chromPosToInt_dict, cellToClusterDict) #building the G' matrix
GprimeDf = GprimeDf.fillna(0)

print("================== G' matrix =========================")
print(GprimeDf)

cell_pos_dict = cell_pos_SNV(initial_matrix, posSet)
filtered_cell_pos_dict = filter_cell_pos_SNV(cell_pos_dict, GprimeDf, clusterToCellDict, posSet)
simulated_input_dict = simulated_input_tree(clusterToCellDict, timeToCellIdDict, filtered_cell_pos_dict)
print(" Simulated input dict =================")
print(simulated_input_dict)
with open(args.output, 'w') as cc:
    json.dump(simulated_input_dict, cc,  indent=4)
