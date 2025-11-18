# Author: Rituparna Khan

import argparse
import numpy as np
import math

''' Get the missing rate for each time point. '''
# Input is cells in each time point and noisy D matrix
def timepoint_missingRate(tpCells, D_matrix):
    tp_MR = {} # Missing rate at each timepoint
    #print(" Timepoint missing rate ")
    noOfMut = len(D_matrix[0])
    for tp in range(len(tpCells)):
        missingEntries = 0
        noOfCells = len(tpCells[tp])
        print(" No. of cells ",noOfCells)
        for i in tpCells[tp]:
            for j in range(noOfMut):
                if D_matrix[int(i)][j] == 3:
                    missingEntries = missingEntries + 1
        missingRate = missingEntries / (noOfCells * noOfMut)
        if missingRate == 0:
            missingRate = float(1e-6)
        tp_MR[tp] = missingRate
    print(" Timepoint missing rate ",tp_MR)
    return tp_MR

# Input is cells in each timepoint file.
# Output a dictionary with timepoint as key and cells as value.
def cellTimepoints(tpCellsFile):
    tcf = open(tpCellsFile, 'r')
    tcf_lines = tcf.readline().rstrip('\n')
    tpCells = []
    while(tcf_lines != ""):
        tp = (tcf_lines.split('\t')[0]).split('_')[1]
        cells = (tcf_lines.split('\t')[1]).split(';')
        tpCells.append(cells)
        tcf_lines = tcf.readline().rstrip('\n')
    #print(tpCells)
    return tpCells

''' Returns the value of P_Dij_gj or the likelihood from Eq. 1. 
From maximum likelihood of the cluster genotype we only consider the 0 or 1 and not the unknowns or 3. '''
# Input is the cell mutation genotype, cluster mutation genotype, alpha, beta and missing rate for each time point.
def prob_cell_j_k(cellMut, clusterMut, alpha, beta, missingRate):
    if clusterMut == 3: # Here 3 is treated as unknowns
        return 0.33
    #print(cellMut," ",clusterMut)
    if (cellMut == 3 and clusterMut == 1) or (cellMut == 3 and clusterMut == 0):
        prob = missingRate
    if cellMut == 1 and clusterMut == 1:
        prob = (1 - missingRate) * (1 - beta)
    if cellMut == 0 and clusterMut == 1:
        prob = (1 - missingRate) * beta
    if cellMut == 1 and clusterMut == 0:
        prob = (1 - missingRate) * alpha
    if cellMut == 0 and clusterMut == 0:
        prob = (1 - missingRate) * (1 - alpha)
    return prob

''' Get the prior value based on value of consensus genotype of subclone k for mutation j. '''
# Input is cluster mutation genotype and missing rate for that time point
def prior_cg_kj(clusterMut, missingRate):
    if clusterMut == 3:
        return missingRate
    if clusterMut == 1 or clusterMut == 0:
        return (1 - missingRate) / 2
    return (1 - missingRate) / 2

''' Get the normalized value or P(D_j). The input will be cells likelihood cells_cj_cgkj,
cluster mutation value, probability of all cells mutation given cluster mutation cgkj_cij
and missingRate to calculate the prior. '''
# Input is cells likelihood cells_cj_cgkj, cluster mutation genotype, probability of all cells mutation given cluster mutation cgkj_cij, missingRate
def normalized_P_Dj(cells_cj_cgkj, clusterMut, cgkj_cij, missingRate):
    # For each mutation we will have the prob_cell_j_k and prior_cg_kj.
    # We just sum up all possible combinations when cg = 0, 1 or unk.
    P_Dj = 0
    possible_cg_kj = [0,1,3]

    for cg_kj in possible_cg_kj:
        if cg_kj == clusterMut: # we use this to not re-calculate this value
            P_cgkj_cij = cgkj_cij
        else:
            P_cgkj_cij =  cells_cj_cgkj * prior_cg_kj(cg_kj, missingRate)
        P_Dj = P_Dj + P_cgkj_cij

    return P_Dj

''' Calculate the probability of the cluster for each timepoint. '''
# Input is cells in each cluster, genotype of the cluster, noisy D matrix, alpha, beta and missing rate for each time point.
def calc_cluster_prob(cluster_cells, cluster_genotype, D_matrix, tp_alpha, tp_beta, tp_missingRate):
    tp_cluster_prob = {}
    total_mutations = len(D_matrix[0])
    print(" Total mutations ",total_mutations)
    for tp, cg in cluster_genotype.items():
        alpha = tp_alpha[tp]
        beta = tp_beta[tp]
        missingRate = tp_missingRate[tp]
        for cluster, gen in cg.items():
            cells = cluster_cells[tp][cluster] # cells in each cluster
            sum_log_P_cgk_ck = 0
            #avg_sum_log_P_cgk_ck = 0
            # j = mutations, i = cell, k = subclones
            #print(" Timepoint ",tp," cluster ",cluster," cells ",cells)
            for j in range(len(gen)):
                if gen[j] == 3: # Skip the missing values
                    continue
                cells_Cj_CGkj = 1 # Cj = all cells of jth mutation, CGkj = consensus genotype of jth mutation of subclone k
                for i in cells:
                    if D_matrix[i][j] == 3: # Skip the missing values
                        continue
                    P_Cij_CGkj = prob_cell_j_k(D_matrix[i][j], gen[j], alpha, beta, missingRate) # likelihood
                    cells_Cj_CGkj = cells_Cj_CGkj * P_Cij_CGkj # multiply the likelihood of the cells

                #print(" Probability of cells ",cells_Cj_CGkj)
                prior_CGkj = prior_cg_kj(gen[j], missingRate)
                #print(" Prior ",prior_CGkj)
                # Numerator of the equation 1
                P_CGkj_Cj = cells_Cj_CGkj * prior_CGkj
                # Normalizer
                P_Dj = normalized_P_Dj(cells_Cj_CGkj, gen[j], P_CGkj_Cj, missingRate)
                #print(" Normalizer ",P_Dj)
                if P_Dj == 0:
                    P_CGk_Ck = 0
                else:
                    P_CGk_Ck = P_CGkj_Cj / P_Dj
                #print(j," mutation prob ",P_CGk_Ck)
                if P_CGkj_Cj == 0:
                    P_CGkj_Cj = float(1e-6)
                with np.errstate(divide='ignore', invalid='ignore'):
                    sum_log_P_cgk_ck = sum_log_P_cgk_ck + np.log(P_CGkj_Cj)
                    #if math.isinf(sum_log_P_cgk_ck):
                    #    print("P_CGkj_Cj value ",P_CGkj_Cj)
                #print(" Inside mutation prob ",sum_log_P_cgk_ck)

            #print(" Prob sum ",sum_log_P_cgk_ck)
            avg_sum_log_P_cgk_ck = (sum_log_P_cgk_ck / total_mutations) / len(cells)
            print("Timepoint ",tp," cluster ",cluster," Cluster probability ",avg_sum_log_P_cgk_ck," # of cells ",len(cells))
            # Save the cluster probability for each timepoint
            tp_cluster_prob[str(tp)+"_"+str(cluster)] = avg_sum_log_P_cgk_ck

    return tp_cluster_prob

# Input is the D matrix file.
# Output is a list of the genotypes with indices corresponding to cells.
def readDMatrix(Dfile):
    dfile = open(Dfile,'r')
    df_line = dfile.readline().rstrip('\n')
    D_matrix = []
    while(df_line != ""):
        mut = df_line.split('\t')
        mut = [int(i) for i in mut]
        D_matrix.append(mut)
        df_line = dfile.readline().rstrip('\n')
    print(" No. of cells ",len(D_matrix)," mut ",len(D_matrix[0]))
    return D_matrix

# Input is the BnpC's log file.
# Output is the FP FN values at each timepoint.
def readFPFNvalues(logF):
    tp_alpha = {}
    tp_beta = {}
    epsilon = float(1e-6)
    print(" Epsilon ",epsilon)
    for i in range(len(logF)):
        lfile = open(logF[i],'r')
        lf_line = lfile.readline().rstrip('\n')
        lf_line = lfile.readline().rstrip('\n')
        lf_arr = lf_line.split('\t')
        fn_rate = lf_arr[3]
        fp_rate = lf_arr[5]
        tp_alpha[i] = float(fp_rate) + epsilon
        tp_beta[i] = float(fn_rate) + epsilon
    print(" TP alpha ",tp_alpha)
    print(" TP Beta ",tp_beta)
    return tp_alpha, tp_beta
