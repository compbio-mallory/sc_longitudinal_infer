import numpy as np
import pandas as pd
from sklearn.metrics.pairwise import nan_euclidean_distances

"""
Function to convert input file/data into binary matrix
parameter: 
in_file : input file
"""
def load_data(in_file, seperator):
    lines = []

    # Get first fine lines to determine if col/row names are provided
    with open(in_file, 'r') as f:
        for i in range(5):
            lines.append(f.readline().strip())
    if seperator == '\t':
        if lines[0].count('\t') > lines[0].count(' '):
            sep = '\t'
        else:
            sep = ' '
    else:
        if lines[0].count(',') > lines[0].count(' '):
            sep = ','
        else:
            sep = ' '

    header_row = False
    header_line = ''
    for el in lines[0].split(sep):
        try:
            el_float = int(el)
        except ValueError:
            if el == ' ':
                continue
            header_row = True
            header_line = lines.pop(0)
            break    
        else:
            if el_float not in [0, 1, 2, 3]:
                header_row = True
                header_line = lines.pop(0)
                break

    index_col = False
    for i, line in enumerate(lines):
        first_el = line.split(sep)[0]
        try:
            first_el_flt = int(first_el)
        except ValueError:
            if first_el == ' ':
                continue
            index_col = True
            break
        else:
            if first_el_flt not in [0, 1, 2, 3]:
                index_col = True
                break

    if index_col and header_row:
        col_types = dict([(j, str) if i == 0 else (j, int) \
            for i,j in enumerate(header_line.split(sep))])
        df = pd.read_csv(in_file, sep=sep, index_col=0, header=0, dtype=col_types)
    elif index_col:
        col_types = dict([(i, str) if i == 0 else (j, int) \
            for i in range(len(lines[0].split(sep)))])
        df = pd.read_csv(in_file, sep=sep, index_col=0, header=None, dtype=col_types)
    elif header_row:
        df = pd.read_csv(in_file, sep=sep, index_col=None, header=0, dtype=int)
    else:
        df = pd.read_csv(in_file, sep=sep, index_col=None, header=None, dtype=int)
    df.replace(3, np.nan, inplace=True)
    # replace homozygos mutations with heterozygos
    df.replace(2, 1, inplace=True)
    return df.values

"""
Function to calucalte Eucledian distance i.e the distance between D and G"
parameter: 
D: observed input matrix
G": true genotype matrix
"""
def calculate_euclidean_distance(G, D):
    return nan_euclidean_distances(G, D)

"""
Function to calucalte E from equation2 in our writeup 
parameter: 
eta_coeff: constant coefficient
lambda_coeff: constant coefficient
sub_clones: number of clusters 
D: observed input matrix
G": true genotype matrix
"""
def calculate_E(eta_coeff, lambda_coeff, sub_clones, G, D):
    distance = calculate_euclidean_distance(G, D)
    # print('The eucliden disatnce between G_" and D is:'.format(distance))
    print('The eucliden disatnce between G" and D is:', distance)
    l1_norm = np.linalg.norm(distance, 1)
    error = (eta_coeff * l1_norm) + (lambda_coeff * sub_clones)
    # print('The error E is:'.format(error))
    print('The error E is:', error)

eta_coeff = 1
lambda_coeff = 0
sub_clones = 1

D = load_data('snv.tsv', '\t')
#D = pd.read_csv('../scg/examples/snv.tsv', sep='\t', index_col=0)
#print(D)
#D = D.replace(3, np.nan)
#D = D.replace(2, 1)
print("D", D)
G = load_data('../Gyanendra_code/GDoublePrimeDf.csv', ',')
print("G", G)

calculate_E(eta_coeff, lambda_coeff, sub_clones, G, D)

# checking the implementation is correct using smaller matrices
A = np.array([[1,2,3],[2,3,4],[0,1,2]])
print("A", A)
B = np.array([[1,2,3],[4,3,2]])
print("B", B)
M = A.shape[0]
N = B.shape[0]

A_dots = (A*A).sum(axis=1).reshape((M,1))*np.ones(shape=(1,N))
print("A_dots", A_dots)

B_dots = (B*B).sum(axis=1)*np.ones(shape=(M,1))
print("B_dots", B_dots)

D_squared =  A_dots + B_dots -2*A.dot(B.T)
print("2*A.dot(B.T)", 2*A.dot(B.T))
print("D_squared", D_squared)

zero_mask = np.less(D_squared, 0.0)
D_squared[zero_mask] = 0.0
distance_1 = np.sqrt(D_squared)
print("Distance between A and B is:", distance_1)
print(" ============== ")
print("L1 norm ", np.linalg.norm(distance_1, 1))
distance_2 = nan_euclidean_distances(A, B)
print("Distance between A and B using sklearn", distance_2)
