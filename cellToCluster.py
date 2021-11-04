"""
For now the output file cluster_posteriors.tsv shoud be in the code directory
"""
import pandas as pd
import json
from collections import defaultdict
  
tsv_data = pd.read_csv('../no_doublet/cluster_posteriors.tsv', sep='\t', index_col=0)

# get the column of max values in every row and make new column with max values
tsv_data['clusters'] = tsv_data.idxmax(axis=1)

# converting cell_id index to column in dataframe
tsv_data['cell_id'] = tsv_data.index
#print(tsv_data)

# grouped_df = tsv_data.groupby("clusters")
# print(grouped_df)
# grouped_lists = grouped_df["cell_id"].apply(list)
# grouped_lists = grouped_lists.reset_index()
# print(grouped_lists)

# as discussion with Ritu we will be storing the result of dataframe as key value-pairs with key
# being the cluster and values being the cell_id
# a = dict(zip(tsv_data.clusters, tsv_data.cell_id))
# print(a)
# a = grouped_lists.set_index('clusters').T.to_dict('list')
# for key, value in a.items():
#     print(key, '->', value)

cellIdToClusterDict = {}

for idx,row in tsv_data.iterrows():
    cellIdToClusterDict[row['cell_id']] = row['clusters']

#for key in cellIdToClusterDict:
#      print(key, '->', cellIdToClusterDict[key])

with open('cellToCluster.json', 'w') as cc:
    json.dump(cellIdToClusterDict, cc,  indent=4)
