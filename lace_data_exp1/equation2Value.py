"""
For now the output file cluster_posteriors.tsv shoud be in the code directory
"""
import pandas as pd
import json
from collections import defaultdict
import argparse
from sklearn.metrics.pairwise import nan_euclidean_distances
import numpy as np
import math

def convertToPandas(fileName):
    pd_fromFile = pd.read_csv(fileName, sep='\t', index_col=0)
    return pd_fromFile

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

''' Function to get the G' matrix. Get it from SCG genotype_posterior and replace SNV 2 with 1.'''
def build_GprimeMatrix(gp_table, cellToClusterDict):
    clusterToCellDict = {} # Building a dictionary with cluster IDs as keys and cell IDs as correspondig values
    #print(" GP table ",gp_table)
    #print(" cellToClusterDict ",cellToClusterDict)
    for cellID in cellToClusterDict:
        clusterID = cellToClusterDict[cellID]
        if clusterID in clusterToCellDict:
            temp_list = clusterToCellDict[clusterID];
            temp_list.append(cellID)
        else:
            clusterToCellDict[clusterID] = [cellID]
    #print("Cluster to cell Dict ",clusterToCellDict)

    GprimeDf = pd.DataFrame()
    for cluster in clusterToCellDict:
        for row in gp_table.itertuples():
            if cluster == str(row.Index) and getattr(row,'event_type') == 'snv' and getattr(row,'probability') > 0.99:
                #prob = getattr(row,'probability')
                #if math.ceil(prob) == 1:
                    pos = getattr(row,'event_id')
                    SNV = getattr(row,'event_value')
                    if SNV == 2.0:
                        SNV = 1.0
                    GprimeDf.loc[cluster,pos] = SNV

    return GprimeDf

''' Function to build the G'' matrix from the output of SCG. '''
def build_GDoublePrimeMatrix(gp_table, cellToClusterDict):
    initial_matrix = pd.DataFrame()
    for cell in cellToClusterDict: # key is the Cell ID so we check for each cell ID its SNV types on a specific position
        clusterID = cellToClusterDict[cell]
        for row in gp_table.itertuples():
            if clusterID == str(row.Index) and getattr(row,'event_type') == 'snv' and getattr(row,'probability') > 0.80:
                #prob = getattr(row,'probability')
                #if math.ceil(prob) == 1:
                    pos = getattr(row,'event_id')
                    SNV = getattr(row,'event_value')
                    initial_matrix.loc[cell,pos] = SNV
        #for index, row in gp_table.iterrows(): 
        #    if clusterID == str(index) and row['event_type'] == 'snv': # Choose the cluster the cell ID belongs to and get the position(event_id) and SNV(event_value) to fill the matrix
        #        pos = row['event_id']
        #        SNV = row['event_value']
                #print("Event value .. ",event_value)
        #        initial_matrix.loc[cell,pos] = SNV
    return initial_matrix

def convertInputMatrixToBinary(D_matrix):
    D_df = pd.read_csv(D_matrix, sep='\t', index_col=0)
    D_df = D_df.replace(2, 1)
    D_df = D_df.replace(3, np.nan)
    return D_df

def findEtaValue(D_df, GDoublePrimeMatrix, cellToClusterDict):
    cellID_metrics= {}
    for cell in cellToClusterDict:
        #print(cell)                                                                                                                                                                                                                                                                                                                                                        
        cell_list = []
        for col in GDoublePrimeMatrix.columns:
            #print(col)                                                                                                                                                                                                                                                                                                                                                     
            if D_df.loc[cell][col] == 1 and GDoublePrimeMatrix.loc[cell][col] == 1.0: #TP = true positive
                cell_list.append("TP")
            elif D_df.loc[cell][col] == 1 and GDoublePrimeMatrix.loc[cell][col] == 0.0: # FP = false positive
                cell_list.append("FP")
            elif D_df.loc[cell][col] == 0 and GDoublePrimeMatrix.loc[cell][col] == 1.0: # FN = false negative
                cell_list.append("FN")
            elif D_df.loc[cell][col] == 0 and GDoublePrimeMatrix.loc[cell][col] == 0.0: # TN = true negative
                cell_list.append("TN")
        cellID_metrics[cell] = cell_list

    FPcount = 0
    FNcount = 0
    TPcount = 0
    TNcount = 0
    for cellID in cellID_metrics:
        metric_list = cellID_metrics[cellID]
        for metric in metric_list:
            if metric == "FP":
                FPcount = FPcount+1
            elif metric == "FN":
                FNcount = FNcount+1
            elif metric == "TP":
                TPcount = TPcount+1
            elif metric == "TN":
                TNcount == TNcount+1
    print("FP count ",FPcount, " FN count ",FNcount, " TP count ",TPcount, " TN count ",TNcount)
    FPrate = FPcount/(FPcount+TPcount)
    FNrate = FNcount/(FNcount+TNcount)
    #print("FP rate ",FPrate, " FN rate ", FNrate)
    eta = 1 - (FPrate + FNrate)/2
    return eta

def calculate_euclidean_distance(G, D):
    return nan_euclidean_distances(G, D)

# G is G'' here
def calculate_E(eta_coeff, lambda_coeff, sub_clones, G, D):
    distance = calculate_euclidean_distance(G, D)
    #print("Distance matrix shape ",distance.shape)
    # print('The eucliden disatnce between G_" and D is:'.format(distance))                                                                                                                                                                                                                                                                                                 
    #print('The eucliden distance between G" and D is:', distance)
    l1_norm = np.linalg.norm(distance, 1)
    error = (eta_coeff * l1_norm) + (float(lambda_coeff) * int(sub_clones))
    # print('The error E is:'.format(error))                                                                                                                                                                                                                                                                                                                                
    #print("The value of Eq. 2", error)
    return error

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (similar to input to SCG)")
parser.add_argument("-gp", "--genotype",dest ="genotype", help="Genotype posterior from SCG")
parser.add_argument("-cp", "--cluster",dest = "cluster", help="Cluster posterior from SCG")
parser.add_argument("-lambda_coeff", "--lambda_coeff",dest = "lambda_coeff", help="Lambda co-efficient of Eq. 2")
parser.add_argument("-subclones", "--subclones",dest = "subclones", help="No of clusters used in SCG")
parser.add_argument("-cellCluster", "--cellCluster",dest = "cellCluster", help="Cells assigned to which clusters with SCG")
args = parser.parse_args()

#gp_path = input("Enter the path to genotype_posteriors.tsv from SCG: ") # Path now is ../no_doublet/genotype_posteriors.tsv
gp_path = args.genotype
gp_table = convertToPandas(gp_path) # gp_table is the dataframe with values from genotype_posteriors.tsv

#cp_path = input("Enter the path to cluster_posteriors.tsv from SCG: ") # Path now is ../no_doublet/cluster_posteriors.tsv
cp_path = args.cluster
cp_table = convertToPandas(cp_path) # cp_table is the dataframe with values from cluster_posteriors.tsv
cp_table['cell_id'] = cp_table.index # cell IDs are the index here

# The dictionary where cell IDs are keys and cluster Ids are the values
with open(args.cellCluster) as json_file:
    cellToClusterDict = json.load(json_file)

GDoublePrimeDf = build_GDoublePrimeMatrix(gp_table, cellToClusterDict) # Building the G'' matrix
#print("=================== G'' MATRIX ========================")
#print(GDoublePrimeDf)
GDoublePrimeDf.to_csv('GDoublePrimeDf.tsv', sep='\t')

GprimeDf = build_GprimeMatrix(gp_table, cellToClusterDict) #building the G' matrix
#print("================== G' matrix =========================")
#print(GprimeDf)

#GDoublePrimeDf = build_GDoublePrimeMatrix(gp_table, GprimeDf, cellToClusterDict) # building the G'' matrix
#print("================== G'' matrix ====================")
#print(GDoublePrimeDf)
# Uncomment the below line to save G'' matrix as a CSV file.
#GDoublePrimeDf.to_csv('GDoublePrimeDf.tsv', sep='\t')

D_matrix_df = convertInputMatrixToBinary(args.input)
eta = findEtaValue(D_matrix_df, GDoublePrimeDf, cellToClusterDict)
print(" Eta value ",eta)
#print(D_matrix_df.shape)

e_value = calculate_E(eta, args.lambda_coeff, args.subclones, GDoublePrimeDf, D_matrix_df)
print(e_value)
