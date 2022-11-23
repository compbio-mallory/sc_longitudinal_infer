import pandas as pd
import argparse
import numpy as np

def read_file(tsvFile,h1):
    if h1 == "false":
        df = pd.read_csv(tsvFile, sep='\t', header=None)
    else:
        df = pd.read_csv(tsvFile, sep='\t', index_col=0)
    return df

#def findTrueClusterNo(GT_df):
#    GT_arr = GT_df.to_numpy()
#    cellToGenotype = {}
#    for x, y in np.ndindex(GT_arr.shape):
        #print(GT_arr[x])
#        cellToGenotype[x] = GT_arr[x]
    
#    sameCellGT_count = {} # cells with same genotype are counted and saved.
#    cellID_list = []
#    for cellID, gt in cellToGenotype.items():
#        count = 1
#        for cellID1, gt1 in cellToGenotype.items():
#            if cellID == cellID1 or cellID in cellID_list:
#                continue
#            if (gt==gt1).all():
#                count = count+1
#                cellID_list.append(cellID1)
#        sameCellGT_count[cellID] = count
#    print(" True cluster size ",len(sameCellGT_count))


def calculate_metric_from_vectors(CG_df,GT_df):
    CG_arr = CG_df.to_numpy()
    GT_arr = GT_df.to_numpy()
    print(CG_arr)
    print(" ============ ")
    print(GT_arr)
    correctNos = 0
    correct1s = 0
    correct0s = 0
    total_1s = 0
    total_0s = 0
    num_cells, num_mutation = CG_arr.shape
    total_entries = num_cells * num_mutation
    for x, y in np.ndindex(CG_arr.shape):
        if CG_arr[x,y] == GT_arr[x,y]:
            correctNos = correctNos+1
        if CG_arr[x,y] == 1 and GT_arr[x,y] == 1:
            correct1s = correct1s+1
        if CG_arr[x,y] == 0 and GT_arr[x,y] == 0:
            correct0s = correct0s+1

    for x, y in np.ndindex(GT_arr.shape):
        if GT_arr[x,y] == 1:
            total_1s = total_1s+1
        if GT_arr[x,y] == 0:
            total_0s = total_0s+1

    print(" All correct entries ",correctNos," total 1s in GT ",total_1s," total 0s in GT ",total_0s)
    print(" All correct 1s ",correct1s)
    print(" All correct 0s ",correct0s)
    print(" All entries ",total_entries)
    accuracy = (correctNos/total_entries)*100
    sensitivity = (correct1s/total_1s)*100
    specificity = (correct0s/total_0s)*100
    print(" Accuracy ",accuracy)
    print(" Sensitivity ",sensitivity)
    print(" Specificity ",specificity)
    TP = correct1s
    FN = total_0s - correct0s
    print(" False negative ",FN)
    recall = TP / (TP + FN)
    print(" Recall ",recall)

parser = argparse.ArgumentParser()
parser.add_argument("-cg", "--cg",dest ="cg", help="Consensus genotype matrix")
parser.add_argument("-gtG", "--gtG",dest ="gtG", help="Ground truth matrix")
parser.add_argument("-h1", "--header",dest ="header", help="Header present or absent")
#parser.add_argument("-doublet", "--doublet",dest ="doublet", help="Doublet data")
#parser.add_argument("-doubletFile", "--doubletFile",dest ="doubletFile", help="Doublet info file")
args = parser.parse_args()

if args.header == "false":
    CG_df = read_file(args.cg,"false")
else:
    CG_df = read_file(args.cg,"true")

GT_df = read_file(args.gtG,"true")

#if args.doublet == "true":
#    doublet_cells_list = get_doublet_cells(args.doubletFile)
#    doublet_GT_CG(GT_df, CG_df, doublet_cells_list)

print(" CG =============================================== ")
print(CG_df)
print(" GT =============================================== ")
print(GT_df)
print(" ==================================================== ")
calculate_metric_from_vectors(CG_df,GT_df)
#findTrueClusterNo(GT_df)
