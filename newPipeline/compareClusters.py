import argparse
import sys
import random

# ================= Step 1 of pipeline ========================
# Input is assignment.txt file and output is dictionary with cell as key and cluster as value.
def get_cell_cluster(assignment):
    with open(assignment,'r') as f:
        lines = f.readlines()

    assignment_list = ((lines[1].split('\t'))[2]).split(" ")
    cell_cluster = {}
    print(len(assignment_list))
    cell_count = 0
    for i in assignment_list:
        if '\n' in i:
            i = i.replace("\n","")
        cell_cluster[cell_count] = i
        cell_count = cell_count+1

    return cell_cluster

# Input is the bnpc assignment.txt and startIndex is the cell number index
# Output is a dictionary with key: cluster and cells: value
def get_tp_cell_cluster(assignment, cells):
    with open(assignment,'r') as f:
        lines = f.readlines()

    assignment_list = ((lines[1].split('\t'))[2]).split(" ")
    cell_cluster = {}
    print(len(assignment_list))
    cells = [int(i) for i in cells]
    #cell_count = startIndex #[i for i in range(28,64)]
    for i in range(len(assignment_list)):
        if '\n' in assignment_list[i]:
            assignment_list[i] = assignment_list[i].replace("\n","")
        #cell_cluster[cell_count] = i
        #cell_count = cell_count+1
        cell_cluster[cells[i]] = assignment_list[i]
    return cell_cluster

# Input is BnpC assignment.txt across all timepoints, startIndex starting cell number for that timepoints,
# endIndex ending cell number for that timepoint.
''' Get the stratified timepoint clusters from the whole assignment. '''
def assignment_onTimepoints(whole_assignment, tp_cells):
    stratified_assignment = {}
    #print(" start index ",startIndex," end index ",endIndex)
    #tp_cells = [i for i in range(startIndex,endIndex+1)]
    tp_cells = [int(i) for i in tp_cells]
    print(tp_cells)
    print(" whole_assignment ",whole_assignment)
    tp_assignment = {key:whole_assignment[key] for key in tp_cells}
    #for cell, cluster in whole_assignment.items():
    #print(" stratified timepoint assignment ",tp_assignment)
    return tp_assignment

''' returns an array of cells belonging to the clusters '''
def groupClusters(cell_cluster):
    cluster_cell = {}
    for cell, cluster in cell_cluster.items():
        if cluster in cluster_cell:
            temp_list = cluster_cell[cluster]
            temp_list.append(cell)
            cluster_cell[cluster] = temp_list
        else:
            cluster_cell[cluster] = [cell]

    print(list(cluster_cell.values()))
    return list(cluster_cell.values())

# Input is list of lists
# Returns a dictionary with index indicating cluster nos
def listToDict(clusterList):
    clusterDict = {}
    for i in range(len(clusterList)):
        clusterDict[i] = clusterList[i]
    return clusterDict

# Input is an array of assignment.txt files from BnpC run for each timepoint, assignment.txt file for BnpC run across all timepoints,
# cells in each timpoint. Output is a dictionary with cells clustered in each timepoint.
def getNewClusters(tp_clustersF, all_clustersF, tpCells):
    #all_clusters = get_cell_cluster(all_clustersF)
    new_tp_clusters = {}
    if len(tp_clustersF) != len(tpCells):
        sys.exit("Your input number of timepoint clusters doesn't match the number of timepoints.")
    #timepoints = [i for i in range(len(tp_clustersF))]
    for i in range(len(tp_clustersF)):
        tp_cluster = get_tp_cell_cluster(tp_clustersF[i], tpCells[i])
        grouped_clusters = groupClusters(tp_cluster)

        #allClusters_tp = assignment_onTimepoints(all_clusters, tpCells[i])
        #grouped_allClusters_tp = groupClusters(allClusters_tp)
        #new_clusters = checkClusterMembership(grouped_clusters, grouped_allClusters_tp)
        #new_tp_clusters[i] = listToDict(new_clusters)
        new_tp_clusters[i] = listToDict(grouped_clusters)

    print(" New clusters ",new_tp_clusters)
    #for tp, clus in new_tp_clusters.items():
    #    print("Timepoint ",tp)
    #    for c, cc in clus.items():
    #        print("Cluster ",c," cells ",len(cc))
    return new_tp_clusters

#[[C0-C3], [C4-C7], [C8-C10]] M2: [[C0-C1], [C2-C3], [C4-C7], [C8-C10]]
# Input is cells from clustering per timepoint and across timepoint.
# Output is new clusters with check of 2 cells belonging to same clusters. 
def checkClusterMembership(tp_cells, all_cells):
    #tp_cells = [[0, 1, 2], [3, 4, 5, 6]]
    #all_cells = [[0, 1, 2], [3, 4]]
    #tp_cells = [[0,1,2,3], [4,5,6,7], [8,9,10]]
    #all_cells = [[0,1], [2,3], [4,5,6,7], [8,9,10]]
    print(" TP clusters ",tp_cells)
    print(" Stratified clusters ",all_cells)
    new_cluster = []
    for cell_i in tp_cells:
        for cell_j in all_cells:
            if cell_i == cell_j:
                new_cluster.append(cell_i)
                #print("Same cluster ",new_cluster)
            else:
            #if cell_i != cell_j:
            #    continue
                temp_cluster = []
                for i in range(len(cell_i)):
                    for j in range(len(cell_j)):
                        if cell_i[i] == cell_j[j]:
                            temp_cluster.append(cell_i[i])
                #print("Temp cluster ",temp_cluster)
                if temp_cluster != []:
                    new_cluster.append(temp_cluster)
            #if (M1(cell_i, cell_j) && M2(cell_i, cell_j))
            #    M3(cell_i, cell_j) 

    print("New clusters after merging ",new_cluster)
    return new_cluster

# ================ End of Step 1 ===============
# ================ Start of Step 2 ================

# Input is cells belonging to a cluster, D_matrix, alpha and beta for specific timepoint
# Output returns the cluster genotype after checking for each mutation in each cluster independently
def calcMaxLikelihood(cells, D_matrix, alpha, beta):
    cluster_genotype = []
    for j in range(len(D_matrix[0])): # j is the mutation
        C_k0 = 1
        C_k1 = 1
        count3 = 0 # count no. of 3s
        count0 = 0 # count no. of 0s
        count1 = 0 # count no. of 1s
        # There should be a condn where all D_ij == 3 then append 3 and not check likelihood
        for c in cells:
            #print(" D_matrix ",D_matrix[c])
            if D_matrix[c][j] == 3:
                count3 = count3+1
                continue
            if D_matrix[c][j] == 0:
                C_k0 = C_k0 * (1 - alpha)
                C_k1 = C_k1 * beta
                count0 = count0+1
                #if j == 7 or j == 8 or j == 10:
                #    print("Inside D=0  C_k0 ",C_k0," C_k1 ",C_k1)
            if D_matrix[c][j] == 1:
                C_k0 = C_k0 * alpha
                C_k1 = C_k1 * (1 - beta)
                count1 = count1+1
                #if j == 7 or j == 8 or j == 10:
                #    print("Inside D=1  C_k0 ",C_k0," C_k1 ",C_k1)
        if j == 13 or j == 18:
            print(j," C_k0 ",C_k0," C_k1 ",C_k1)
        # If there is a tie in no. of 0s and 1s and 3s then it is an unknown
        #if count0 == count1 or count3 == len(cells):
        #    if j == 18:
        #        print("Under unknown condn")
        #    cluster_genotype.append(3) # we input 3 for unknowns as well
        if C_k0 > C_k1:
            cluster_genotype.append(0)
        elif C_k1 >= C_k0:
            cluster_genotype.append(1)
        #else:
        #    if j == 18:
        #        print("Under C_k0 == C_k1")
            # If C_k0 == C_k1 then randomly select 0 or 1
        #    random_gen = random.randint(0, 1)
        #    cluster_genotype.append(random_gen)
    #print(" Cluster genotype ",cluster_genotype)
    return cluster_genotype

# Input is D_matrix and the clusters in each timepoint.
# Output is the updated genotype for each cluster in each timepoint.
def updateTpClusterGen(D_matrix, tp_clusters, tp_alpha, tp_beta):
    tp_cluster_gen = {}
    for tp, clusters in tp_clusters.items():
        alpha = tp_alpha[tp]
        beta = tp_beta[tp]
        cluster_gen = {}
        for cluster, cells in clusters.items():
            print(" Timepoint ",tp," Cluster ",cluster)
            # Print the D_matrix for the cells
            for c in cells:
                print("Cell ",c," ",D_matrix[c])

            if len(cells) == 1:
                cluster_gen[cluster] = D_matrix[cells[0]]
                print(cluster_gen[cluster])
            else:
                cluster_gen[cluster] = calcMaxLikelihood(cells, D_matrix, alpha, beta)
                print(cluster_gen[cluster])

        #print("Timepoint ",tp," Cluster ",cluster," gen ",cluster_gen)
        tp_cluster_gen[tp] = cluster_gen
    return tp_cluster_gen

# ================ End of Step 2

# ================ Start of Step 3
# Input is two cluster genotypes
# Returns False if there is a mismatch in terms of 1 or 0. Otherwise True for a match.
def checkSim(c1,c2):
  mismatches = 0
  for i in range(len(c1)):
      if c1[i] != c2[i] and (c1[i] != 3 and c2[i] != 3):
        #print(" Not match ")
        mismatches = mismatches+1
  if mismatches > 0:
    return False
  else:
    return True

# Input is a cluster and list of merged clusters
# Outputs True if a cluster is already merged
def checkIfClusterMerged(cluster, finalClusters):
  if finalClusters == []:
    return False
  for fc in finalClusters:
    if cluster in fc:
      return True

# Input is two cluster genotypes
# Output is the updated cluster genotype comparing the two genotypes.
def updateGenotype(gen1,gen2):
  finalGen = []
  for i in range(len(gen1)):
    if gen1[i] == gen2[i]:
      finalGen.append(gen1[i])
    elif gen1[i] == 3 and gen2[i] != 3:
      finalGen.append(gen2[i])
    elif gen2[i] == 3 and gen1[i] != 3:
      finalGen.append(gen1[i])
    else:
      finalGen.append(0)
  return finalGen

# Input is dictionary of merged clusters, their updated genotypes and all clusters in a timepoint
# For singleton clusters add it to the updatedCG
def allClusters_updateCG(mergedCluster, updatedCG, clusters):
    for c, cc in mergedCluster.items():
        if c not in updatedCG:
            updatedCG[c] = clusters[c]
    return updatedCG

# Input is the clusters underlying genotype obtained using maximum likelihood in each timepoint
# Output is updated cluster genotype where we merge clusters based on their genotype
def checkClusterGenotype(tp_clustersGen):
    tp_mergedClusters = {}
    tp_updatedCG = {}
    for tp, clustersGen in tp_clustersGen.items():
        finalMergedClusters = []
        mergedCluster = {}
        updatedCG = {}
        for c1, gen1 in clustersGen.items():
            clusterList = [c1]
            if checkIfClusterMerged(c1, finalMergedClusters):
                continue
            for c2, gen2 in clustersGen.items():
                if c1 == c2 or checkIfClusterMerged(c2, finalMergedClusters):
                    continue
                if checkSim(gen1, gen2):
                    if c1 in updatedCG:
                        updatedCG[c1] = updateGenotype(updatedCG[c1], gen2)
                    else:
                        updatedCG[c1] = updateGenotype(gen1, gen2)
                    clusterList.append(c2)

            finalMergedClusters.append(clusterList)
            mergedCluster[c1] = clusterList

        updatedCG = allClusters_updateCG(mergedCluster, updatedCG, clustersGen)
        #print("Timepoint ",tp,"Merged clusters ",mergedCluster)
        print("Timepoint ",tp,"Updated cluster genotype ",updatedCG)
        tp_mergedClusters[tp] = mergedCluster
        tp_updatedCG[tp] = updatedCG
    return tp_mergedClusters, tp_updatedCG

# Input is merged clusters in each timepoint and cells for each clusters.
# Output is reassigned cells from other clusters into the merged cluster.
def reassignCells(tp_mergedClusters, tp_cluster_cells):
    tp_reassignedCells = {}
    for tp, mc in tp_mergedClusters.items():
        cluster_cells = {}
        for c1, c2 in mc.items():
            if len(c2) == 1 and c1 == c2[0]: # If no clusters are merged then no need to reassign cells
                cluster_cells[c1] = tp_cluster_cells[tp][c1]
                continue
            cells = tp_cluster_cells[tp][c1]
            for c in c2:
                cells.extend(tp_cluster_cells[tp][c])
            cluster_cells[c1] = list(set(cells))

        tp_reassignedCells[tp] = cluster_cells
    print(" After reassigning cells for timepoint clusters ",tp_reassignedCells)
    return tp_reassignedCells

# ========================== End of Step 3 ===========================

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
    
#parser = argparse.ArgumentParser()
#parser.add_argument("-tp", "--tp",dest = "tp", help="Timepoint cells genotype.")
#parser.add_argument("-cg", "--cg",dest = "cg", help="Consensus genotype from BnpC.")
# Get the timepoint clusters input as an array because we won't know how many timepoints we will have.
# Get the cells assigned in each timepoint instead of using index.
#parser.add_argument("-tc", "--tclusters", nargs="*", help="Timepoint cluster assignment files.")
#parser.add_argument("-cells","--cells", help="Cells assigned at each timepoint.")
#parser.add_argument("-cluster", "--cluster",dest="cluster", help="Assignment.txt across all timepoints.")
#parser.add_argument("-D","--D", help="D matrix.")
#parser.add_argument("-logFile","--logFile", nargs="*", help="BnpC error file to read FP and FN rates for each timepoint.")
# User should be able to enter FP and FN rates if known while running BnpC and enter those values here.
#args = parser.parse_args()

#tpCells = cellTimepoints(args.cells)
#D_matrix = readDMatrix(args.D)
#tp_alpha, tp_beta = readFPFNvalues(args.logFile)

#tp_cluster_cells = getNewClusters(args.tclusters, args.cluster, tpCells)
#tp_cluster_gen = updateTpClusterGen(D_matrix, tp_cluster_cells, tp_alpha, tp_beta)
#print(" TP cluster gen ",tp_cluster_gen)
#tp_mergedClusters, tp_updatedCG = checkClusterGenotype(tp_cluster_gen)
#print(" TP merged clusters ",tp_mergedClusters)
#print(" TP updated CG ",tp_updatedCG)
#print(" TP cluster cells ",tp_cluster_cells)
#tp_reassignedCells = reassignCells(tp_mergedClusters, tp_cluster_cells)
