import json
import numpy as np
import pandas as pd
import argparse

def find_cols_index(D_matrix):
    D_df = pd.read_csv(D_matrix, sep = '\t')
    cells = list(D_df['cell_id'])
    cols = list(D_df.columns)
    #print(cells)
    #print(cols[1:])
    return cells, cols[1:]

def find_Bmatrix(B_list):
    B_matrix = np.reshape(B_list,(13,13),order='C')
    return B_matrix

def build_GDoubleprimeMatrix(B_matrix, C_dict, cells, chrom_pos):
    initial_matrix = pd.DataFrame(index = cells, columns = chrom_pos)
    cell_to_clones = {}
    t_cell_clones = C_dict['Experiment_1']
    x_cell_clones = C_dict['Experiment_2']
    t_index = 0
    x_index = 0
    for cell in cells: # need to create a dict with cells and their corresponding clones.
        if 'T' in cell:
            cell_to_clones[cell] = t_cell_clones[t_index]
            t_index = t_index+1
        else:
            cell_to_clones[cell] = x_cell_clones[x_index]
            x_index = x_index+1   

    #print(cell_to_clones)
    for cell, clone_no in cell_to_clones.items():
        #print(cell," ",clone_no)
        B_matrix_list = B_matrix[clone_no]
        #print(B_matrix_list[1:])
        initial_matrix.loc[cell] = B_matrix_list[1:]
    
    #print(initial_matrix)
    return initial_matrix

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (similar to input to SCG)")
args = parser.parse_args()
cells,chrom_pos = find_cols_index(args.input)

B_list = []
C_dict = {}
clones_summary = {}
with open('/gpfs/research/fangroup/rk18g/longitudinal/LACE-UTILITIES/real_data/Breast_Cancer/p494/results/inference.json') as fp:
    lace_dict = json.load(fp)
    B_list = lace_dict['B']
    C_dict = lace_dict['C']
    clones_summary = lace_dict['clones_summary']  
    #print(lace_dict)

print("B ----- ",B_list)
print("C ----- ",C_dict)
print("Clones summary ",clones_summary)
print(" ",cells)
print(chrom_pos)
B_matrix = find_Bmatrix(B_list) #Ignore the last column which will be clone 13 as it doesn't exist.
print(B_matrix)
gDoublePrimeMatrix = build_GDoubleprimeMatrix(B_matrix, C_dict, cells, chrom_pos)
gDoublePrimeMatrix.to_csv('LACE_G_matrix.tsv',sep='\t')
