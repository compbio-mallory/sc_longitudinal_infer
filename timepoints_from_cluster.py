"""
This file to simulate time points for each cells.
"""
import pandas as pd
import json
from collections import defaultdict
  
tsv_data = pd.read_csv('../no_doublet/cluster_posteriors.tsv', sep='\t', index_col=0)
#tsv_data = pd.read_csv('example_cellToCluster.json', sep='\t', index_col=0)

# get the column of max values in every row and make new column with max values
tsv_data['time_point'] = tsv_data.idxmax(axis=1)

# converting cell_id index to column in dataframe
tsv_data['cell_id'] = tsv_data.index
#print(tsv_data)

cellIdToTimePointDict = {}

for idx,row in tsv_data.iterrows():
    time_point = int(row['time_point']) + 1
    if time_point in cellIdToTimePointDict.keys():
        cell_list = cellIdToTimePointDict[time_point]
        cell_list.append(row['cell_id'])
        cellIdToTimePointDict[time_point] = cell_list
    else:
        cellIdList = []
        cellIdList.append(row['cell_id'])
        cellIdToTimePointDict[time_point] = cellIdList

#for key in cellIdToClusterDict:
#      print(key, '->', cellIdToClusterDict[key])

# need to reverse this to include timepoints --> cells
with open('cell_timepoints.json', 'w') as cc:
    json.dump(cellIdToTimePointDict, cc,  indent=4)
