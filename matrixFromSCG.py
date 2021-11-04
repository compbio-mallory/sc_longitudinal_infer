"""
For now the output file cluster_posteriors.tsv shoud be in the code directory
"""
import pandas as pd
import json
from collections import defaultdict

def convertToPandas(fileName):
    pd_fromFile = pd.read_csv(fileName, sep='\t', index_col=0)
    return pd_fromFile

''' To initialize the matrix with 
cell IDs as index and positions as columns'''
def init_rowsColumns(gp_table, cp_table):
    pos_list = gp_table['event_id'].values
    df = pd.DataFrame(index = cp_table.index,columns = pos_list)
    #df = pd.DataFrame(index = range(noOfCells),columns = pos_list)
    return df

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
def build_GprimeMatrix(gp_table, final_matrix, cellToClusterDict):
    clusterToCellDict = {} # Building a dictionary with cluster IDs as keys and cell IDs as correspondig values
    for key in cellToClusterDict:
        clusterID = cellToClusterDict[key]
        if clusterID in clusterToCellDict:
            temp_list = clusterToCellDict[clusterID];
            temp_list.append(key)
        else:
            clusterToCellDict[clusterID] = [key]
    #print(clusterToCellDict.keys())
    posSet = set(gp_table['event_id'])

    GprimeDf = pd.DataFrame()
    # Below code builds the matrix with cluster IDs as indices and positions as columns.
    # Fills it up with voted SNV values after creating a list with the SNV values for each cell ID.
    for cluster in clusterToCellDict:
        cell_SNV_list = []
        for pos in posSet:
            for cell in clusterToCellDict[cluster]:
                temp_value = final_matrix.loc[cell][pos].to_string()
                #print(temp_value.split()[1])
                cell_SNV_list.append(temp_value.split()[1])
            #Include voting results here
            voted_SNV = voting_algo(cell_SNV_list)
            GprimeDf.loc[cluster,pos] = voted_SNV 
            #print("SNV list ::",cell_SNV_list)

    return GprimeDf

''' Function to get the G'' matrix '''
def build_GDoublePrimeMatrix(gp_table, GprimeDf, cellToClusterDict):
    GDoublePrimeDf = pd.DataFrame()
    posSet = set(gp_table['event_id'])
    for cell in cellToClusterDict:
        clusterID = cellToClusterDict[cell]
        for cluster, row in GprimeDf.iterrows():
            if clusterID == str(cluster):
                for pos in posSet:
                    SNV_value = row[pos]
                    GDoublePrimeDf.loc[cell,pos] = SNV_value
    return GDoublePrimeDf

''' Function to build the initial matrix from the output of SCG. This will look similar to matrix D. '''
def final_initial_matrix(gp_table, cellToClusterDict, initial_matrix):
    for cell in cellToClusterDict: # key is the Cell ID so we check for each cell ID its SNV types on a specific position
        clusterID = cellToClusterDict[cell]
        for index, row in gp_table.iterrows(): 
            if clusterID == str(index): # Choose the cluster the cell ID belongs to and get the position(event_id) and SNV(event_value) to fill the matrix
                pos = row['event_id']
                SNV = row['event_value']
                #print("Event value .. ",event_value)
                initial_matrix.loc[key,pos] = SNV
    return initial_matrix
    
gp_table = convertToPandas('../no_doublet/genotype_posteriors.tsv') # gp_table is the dataframe with values from genotype_posteriors.tsv
cp_table = convertToPandas('../no_doublet/cluster_posteriors.tsv') # cp_table is the dataframe with values from cluster_posteriors.tsv
cp_table['cell_id'] = cp_table.index # cell IDs are the index here

# The dictionary where cell IDs are keys and cluster Ids are the values
with open('cellToCluster.json') as json_file:
    cellToClusterDict = json.load(json_file)
    
initial_matrix = init_rowsColumns(gp_table, cp_table) # Get the initial matrix similar to D with values from genotype_posteriors.tsv and cluster_posteriors.tsv (output from SCG)
print("=================== INITIAL MATRIX ========================")
print(initial_matrix)

final_matrix=final_initial_matrix(gp_table, cellToClusterDict, initial_matrix) # Populated matrix similar to D
print("=================== FINAL MATRIX ========================")
print(final_matrix)

GprimeDf = build_GprimeMatrix(gp_table, final_matrix, cellToClusterDict) #building the G' matrix
print("================== G' matrix =========================")
print(GprimeDf)

GDoublePrimeDf = build_GDoublePrimeMatrix(gp_table, GprimeDf, cellToClusterDict) # building the G'' matrix
print("================== G'' matrix ====================")
print(GDoublePrimeDf)
# Uncomment the below line to save G'' matrix as a CSV file.
#GDoublePrimeDf.to_csv('GDoublePrimeDf.csv')
