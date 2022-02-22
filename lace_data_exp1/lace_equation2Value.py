import pandas as pd
import json
from collections import defaultdict
import argparse
from sklearn.metrics.pairwise import nan_euclidean_distances
import numpy as np
import math

def convertInputMatrixToBinary(D_matrix):
    D_df = pd.read_csv(D_matrix, sep='\t', index_col=0)
    D_df = D_df.replace(2, 1)
    D_df = D_df.replace(3, np.nan)
    D_df = D_df.fillna(0)
    return D_df

def findEtaValue(D_df, GDoublePrimeMatrix):
    cellID_metrics= {}
    cells_list = D_df.index.values.tolist()
    for cell in cells_list:
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
parser.add_argument("-inputD", "--inputD",dest ="inputD", help="Input matrix D")
parser.add_argument("-inputG", "--inputG",dest ="inputG", help="Input matrix G")
args = parser.parse_args()

D_matrix_df = convertInputMatrixToBinary(args.inputD)
G_matrix_df = convertInputMatrixToBinary(args.inputG)
print(D_matrix_df)
print(G_matrix_df)
eta = findEtaValue(D_matrix_df, G_matrix_df)
print(eta)
e_value = calculate_E(eta, 0.2, 12, G_matrix_df, D_matrix_df)
print(e_value)

