import numpy as np

"""
Function to calucalte Eucledian distance i.e the distance between D and G"
parameter: 
D: observed input matrix
G": true genotype matrix
"""
def calculate_euclidean_distance(G, D):
    p1 = np.sum(G**2, axis=1)[:, np.newaxis]
    p2 = np.sum(D**2, axis=1)
    p3 = -2 * np.dot(G, D.T)
    return np.round(np.sqrt(p1+p2+p3), 2)

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
D = np.array([[0,1],[0,1],[1,0],[1,0]])
G = np.array([[0,0],[1,1],[1,0],[1,0]])
calculate_E(eta_coeff, lambda_coeff, sub_clones, G, D)