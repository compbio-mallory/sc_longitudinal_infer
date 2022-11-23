"""
For now the output file cluster_posteriors.tsv shoud be in the code directory
"""
import pandas as pd
import json
from collections import defaultdict

def cell_cluster(cluster_posterior,jsonOutputFile):
    tsv_data = pd.read_csv(cluster_posterior, compression='gzip', sep='\t', index_col=0)
    # get the column of max values in every row and make new column with max values
    tsv_data['clusters'] = tsv_data.idxmax(axis=1)
    tsv_data['cell_id'] = tsv_data.index # converting cell_id index to column in dataframe
    cellIdToClusterDict = {}
    for idx,row in tsv_data.iterrows():
        cellIdToClusterDict[row['cell_id']] = row['clusters']

    with open(jsonOutputFile, 'w') as cc:
        json.dump(cellIdToClusterDict, cc,  indent=4)

#cell_cluster('../no_doublet/cluster_posteriors.tsv','cellToCluster.json')
