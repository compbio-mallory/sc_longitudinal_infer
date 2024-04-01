import argparse
import numpy as np
import pydot
import random
import os
import copy
from compareClusters import getNewClusters, updateTpClusterGen, checkClusterGenotype, reassignCells
from calculateClusterProb import calc_cluster_prob, timepoint_missingRate
from treeAlgorithm import Node, buildTree, printTree, connectSameNodes, find_unobservedSubclone, connect_remainingNodes, getEdges
from checkParallelAndBackMut import getParallelMut, getBackMutCount, finalParallelMutEdges, finalBackMutEdges, correctParallelMut, correctBackMut
from visualizeTree import drawTree

''' In a given timepoint calculate the likelihood of the cells for the clusters. '''
# Input is noisy D matrix, alpha beta values for each timepoint, cells to be reassigned, clusters and their genotype in each timepoint
def cell_cluster_likelihood(D_matrix, alpha, beta, cells, cluster_genotype, ignore_cluster, cluster_cells):
    # If there are no more subclones left to be merged in a timepoint then don't check any further
    if len(cluster_genotype) == 1 and ignore_cluster in cluster_genotype:
        return "Done", "Done"

    del(cluster_cells[ignore_cluster]) # Remove the spurious cluster
    del(cluster_genotype[ignore_cluster])
    
    for i in cells:
        p_Di_Ck = 0 # i = cell, j = mutation, k = subclone
        total_mutations = len(D_matrix[0])
        cluster_prob = {}
        for cluster, gen in cluster_genotype.items():
            if cluster == ignore_cluster:
                continue
            ignore_mut = 0
            for j in range(len(gen)):
                if gen[j] == 3 or D_matrix[i][j] == 3:
                    ignore_mut = ignore_mut+1
                    continue
                if D_matrix[i][j] == 0 and gen[j] == 1:
                    L_j = beta
                if D_matrix[i][j] == 1 and gen[j] == 1:
                    L_j = 1 - beta
                if D_matrix[i][j] == 1 and gen[j] == 0:
                    L_j = alpha
                if D_matrix[i][j] == 0 and gen[j] == 0:
                    L_j = 1 - alpha
                p_Di_Ck = p_Di_Ck + np.log(L_j)
            if (total_mutations - ignore_mut) == 0:
                avg_prob = 0
            else:
                avg_prob = p_Di_Ck / (total_mutations - ignore_mut)
            cluster_prob[cluster] = avg_prob
        # select the cluster with max likelihood
        print(" TP cluster ",cluster_prob," cell ",i)
        max_cluster = max(cluster_prob,key=cluster_prob.get)
        print(" cell ",i," max cluster ",max_cluster)
        # Reassigning the cells to the selected cluster
        cluster_cells[max_cluster].append(i)

    return cluster_cells, cluster_genotype

''' After removing spurious subclones, reassigning cells we need to update cluster genotype. '''
def correctClustering(D_matrix, tp_cluster_cells, tp_alpha, tp_beta):
    tp_cluster_gen = updateTpClusterGen(D_matrix, tp_cluster_cells, tp_alpha, tp_beta)
    print(" TP cluster gen ",tp_cluster_gen)
    tp_mergedClusters, tp_updatedCG = checkClusterGenotype(tp_cluster_gen)
    print(" TP merged clusters ",tp_mergedClusters)
    print(" TP updated CG ",tp_updatedCG)
    print(" TP cluster cells ",tp_cluster_cells)
    tp_reassignedCells = reassignCells(tp_mergedClusters, tp_cluster_cells)
    return tp_reassignedCells, tp_updatedCG

''' Infer the tree using timepoint clusters, cells and genotype. '''
# Input is cells in each cluster and cluster genotype
def getTree(tp_cluster_cells, tp_cluster_genotype):
    Tree, clone_node, tp_nodes = buildTree(tp_cluster_cells, tp_cluster_genotype)
    #print("Clone nodes ",clone_node)
    Tree, connectedNodes = connectSameNodes(Tree)
    Tree, tp_usc_nodes = find_unobservedSubclone(Tree, tp_nodes, connectedNodes)
    print("Unobserved subclones in each timepoint ",tp_usc_nodes)
    Tree = connect_remainingNodes(Tree, tp_nodes, tp_usc_nodes)
    Tree = getEdges(Tree)
    printTree(Tree)
    return Tree, clone_node

''' Get details from the Tree and generate plots. '''
def plotTree(Tree, plMut_edges, backMut_edges, opFile, name, sample):
    tp_nodes = {}
    id_cells = {}
    id_newMut = {}
    id_plMut = {}
    id_backMut = {}
    id_edges = {}
    usc_nodes = []
    skip_nodes = []
    for i in range(len(Tree)):
        if Tree[i].pID == -2: # Don't consider this node for plotting
            skip_nodes.append(i)
            continue
        if Tree[i].type == "unobserved":
            usc_nodes.append(Tree[i].id)
        # Add the nodes in each timepoint
        if Tree[i].timepoint in tp_nodes:
            tp_nodes[Tree[i].timepoint].append(i)
        else:
            tp_nodes[Tree[i].timepoint] = [i]
        id_cells[i] = Tree[i].cells
        id_edges[i] = Tree[i].edges
        # Get the new and back mutation for this node
        p_mut = Tree[Tree[i].pID].mutations
        new_mut = set(Tree[i].mutations) - set(p_mut)
        id_newMut[i] = list(new_mut)

    for m, edges in backMut_edges.items():
        for e in edges:
            nodeId = int(e.split('_')[1])
            if nodeId in skip_nodes: # skip adding nodes already removed
                continue
            if nodeId in id_backMut:
                id_backMut[nodeId].append(m)
            else:
                id_backMut[nodeId] = [m]

    #if name == "Final Tree":
    #    print("Back mutation ",id_backMut)

    for p, edges in plMut_edges.items():
        for e in edges:
            if e == '':
                continue
            nodeId = int(e.split('_')[1])
            if nodeId in id_plMut:
                id_plMut[nodeId].append(p)
            else:
                id_plMut[nodeId] = [p]

    plottedTree = drawTree(tp_nodes, id_cells, id_newMut, id_plMut, id_backMut, id_edges, usc_nodes, name, sample)
    #print(plottedTree)
    plottedTree.write_png(opFile)
    #print(" Timepoint nodes ",tp_nodes)
    #print(" ID cells ",id_cells)
    #print(" ID new mut ",id_newMut)
    #print(" ID edges ",id_edges)
    #print(" ID back mut ",id_backMut)

''' Given the probability of all the clones sum up their probabilities. '''
def calculate_treeProb(tp_cluster_prob):
    tree_prob = 0
    for c,p in tp_cluster_prob.items():
        tree_prob = tree_prob + p
    return tree_prob

''' Select only those subclones < upper bound for each timepoint and sort them increasingly with probabilities. '''
# Input is the cluster prob threshold for each timepoint and prob of clusters in each timepoint
def getSpuriousSubclones(tp_thresholds, cluster_prob):
    selected_clusters = {}
    for cluster, prob in cluster_prob.items():
        tp = int(cluster.split("_")[0])
        if prob < tp_thresholds[tp]:
            selected_clusters[cluster] = prob

    # Then sort the selected clusters with increasing probability values
    sort_selected_clusters = dict(sorted(selected_clusters.items(), key=lambda x:x[1]))
    #sort_selected_clusters = dict(sort_selected_clusters)
    #sort_selected_clusters = list(sort_selected_clusters.keys())
    print(" Subclones < upper bound and sorted by increasing probabilities.")
    #print(" Cluster prob ",cluster_prob)
    print(" Sort selected clusters ",sort_selected_clusters)
    return sort_selected_clusters

''' Get weighted thresholds for all timepoints. '''
# Input is prob of clusters in each timepoint and cells in each timepoint
def clusterUpperBound(tp_cluster_prob, tpCells, tp_cluster_cells):
    tp_thresholds = {}
    tp_prob = {}
    
    for tp_cluster, prob in tp_cluster_prob.items():
        tp = int(tp_cluster.split('_')[0])
        cluster = int(tp_cluster.split('_')[1])
        cells = tp_cluster_cells[tp][cluster]
        cluster_prob = prob * len(cells)

        if tp in tp_prob:
            tp_prob[tp] = tp_prob[tp] + cluster_prob
        else:
            tp_prob[tp] = cluster_prob

    for i in range(len(tpCells)):
        #print(" Timepoint ",i," cells ",len(tpCells[i]))
        tp_thresholds[i] = tp_prob[i] / len(tpCells[i])

    print(" Timepoint thresholds ",tp_thresholds)
    # Find the max weighted avg of all time points which will be the upper bound. Return the upper bound.
    #upper_bound = max(tp_thresholds.values())
    #print(" Upper bound of the threshold ",upper_bound)
    return tp_thresholds

''' After correction some subclones might need merging for having same set of mutations. '''
def mergeClones(Tree):
    # Get the time point nodes. For each node in same timepoint check and merge with cells reassigned.
    tp_nodes = {}
    for i in range(len(Tree)):
        if Tree[i].type == "unobserved":
            continue
        tp = Tree[i].timepoint
        if tp in tp_nodes:
            nodes = tp_nodes[tp]
            for n in nodes:
                if set(Tree[i].mutations) == set(Tree[n].mutations):
                    Tree[n].children.extend(Tree[i].children)
                    Tree[n].cells.extend(Tree[i].cells)
                    Tree[Tree[i].pID].children.remove(i)
                    Tree[i].pID = -2
                    Tree[i].children = []
                    Tree[i].edges = ""
        else:
            tp_nodes[tp] = [i]
    print("Tree after merging ===================")
    printTree(Tree)
    return Tree

''' Re-check and correct the Tree after fixing parallel and back mutations. '''
def recheckTree(Tree):
    #printTree(Tree)
    # Check if at any time point two subclones have same set of mutations and merge them.
    #Tree = mergeClones(Tree)
    # Then check after correction if any unobserved subclones are unnecessary
    #nodesToRemove = set()
    for i in range(len(Tree)):
        #print(" Tree i ",i)
        if Tree[i].pID == -2: # Don't check for already removed clones
            continue
        if Tree[i].type == "unobserved":
            if Tree[i].pID == 0:
                if len(Tree[i].children) == 1: # If there is just 1 children then remove this unobserved subclone
                    Tree[Tree[i].children[0]].pID = 0
                    Tree[Tree[i].children[0]].edges = str(Tree[i].pID)+"_"+str(Tree[i].children[0])
                    Tree[Tree[i].pID].children.remove(i)
                    Tree[Tree[i].pID].children.extend(Tree[i].children)
                    Tree[i].pID = -2
                    Tree[i].children = []
                    Tree[i].edges = ""
            elif set(Tree[Tree[i].pID].mutations) == set(Tree[i].mutations):
                #print(Tree[Tree[i].pID].mutations," ",Tree[i].mutations)
                Tree[Tree[i].pID].children.remove(i)
                Tree[Tree[i].pID].children.extend(Tree[i].children)
                for c in Tree[i].children:
                    Tree[c].pID = Tree[i].pID
                    Tree[c].edges = str(Tree[i].pID)+"_"+str(c)
                print(" Nodes to remove ",i)
                Tree[i].pID = -2
                Tree[i].children = []
                Tree[i].edges = ""
                #nodesToRemove.add(i)

    #Treelen = len(Tree)
    #for n in nodesToRemove:
    #    if n >= Treelen:
    #        Tree.remove(Tree[Treelen - 1])
    #    else:
    #        Tree.remove(Tree[n])
    #        Treelen = len(Tree)
    return Tree

''' Recheck final tree before plotting. '''
def recheckFinalTree(Tree):
    # Merge clones having similar set of mutations.
    # Then check for unobserved subclones.
    Tree = mergeClones(Tree)
    Tree = recheckTree(Tree)
    return Tree

''' Update the genotype of the clusters from the corrected Tree. '''
# Input is Tree, cluster genotype in each timepoint, clusters to nodes in Tree.
def genotypeFromTree(Tree, tp_cluster_gen, clone_node, noOfMutations):
    print("Clone node ",clone_node)
    print(" tp_cluster_gen ",tp_cluster_gen)
    for c,n in clone_node.items():
        tp = int(c.split('_')[0])
        cluster = int(c.split('_')[1])
        mut = Tree[n].mutations
        if mut == []:
            genotype = [0] * noOfMutations
        else:
            genotype = [0] * noOfMutations
            for n in mut:
                genotype[n] = 1
        tp_cluster_gen[tp][cluster] = genotype
    print("Genotype ",tp_cluster_gen)
    return tp_cluster_gen

#''' If any parallel edge appear as back mutation then keep that as back mutation. '''
# Input is parallel and mutation egdes
#def filterParallelEdges(parallelEdges, backEdges):
#    for bm in backEdges:
#        if bm in parallelEdges:
#            del(parallelEdges[bm])
#    return parallelEdges

''' Prints the mutation from genotype '''
def printMutFromGenotype(genotype):
    for c,g in genotype.items():
        mut = []
        for i in range(len(g)):
            if g[i] == 1:
                mut.append(i)
        print(str(c)," ",mut)

''' Note the diff between the old and new FP FN rates. '''
def diffFPFN(old_tp_FP, old_tp_FN, new_tp_FP, new_tp_FN):
    tp_diffFP = {}
    tp_diffFN = {}
    for tp in old_tp_FP:
        print(" Timepoint ",tp," Old FP ",old_tp_FP[tp]," Old FN ",old_tp_FN[tp])
        print(" Timepoint ",tp," New FP ",new_tp_FP[tp]," New FN ",new_tp_FN[tp])
        diff_FP = new_tp_FP[tp] - old_tp_FP[tp]
        diff_FN = new_tp_FN[tp] - old_tp_FN[tp]
        print(" Timepoint ",tp," Diff FP ",diff_FP," Diff FN ",diff_FN)
        tp_diffFP[tp] = diff_FP
        tp_diffFN[tp] = diff_FN
    return tp_diffFP, tp_diffFN

''' Check diff between FP FN rates. If it is too high then we don't drop this spurious subclone. '''
def checkDiff_FPFN(tp_diffFP, tp_diffFN, tp, userRate):
    #for tp in tp_diffFP:
    FP = tp_diffFP[tp]
    FN = tp_diffFN[tp]
    if abs(FP) >= userRate or abs(FN) >= userRate:
        return True

''' Select the tree with the highest probability. '''
# Input is D_matrix, alpha and beta in each timepoint, cluster probabilities in each timepoint, cells in each cluster, cluster genotypes
def selectOptimalTree(D_matrix, tp_alpha, tp_beta, tp_MR, tp_cluster_prob, tp_cluster_cells, tp_cluster_genotype, tpCells, k, userRate, plotOp, bnpcRun, sample):
    # Get the current tree_prob
    # sort the cluster_prob from least to highest
    # Loop over the cluster_prob and eliminate the least cluster. Reassign its cells to another cluster based on maximum likelihood. Get the cluster genotype after reassigning cells.
    # Infer the tree with new clustering results
    # Run parallel and back mutation correction algorithm
    # Calculate new tree_prob 
    # If new tree_prob < tree_prob then break. And resulting tree is one with the highest probability.
    Tree, clone_node = getTree(tp_cluster_cells, tp_cluster_genotype)
    #print(" Clone node ",clone_node)
    tree_prob = calculate_treeProb(tp_cluster_prob) 
    # Get the prob threshold for each timepoint
    upper_bound = clusterUpperBound(tp_cluster_prob, tpCells, tp_cluster_cells)
    print("Upper bound ",upper_bound)
    print("Tree prob ",tree_prob)
    # Only subclones to eliminate
    spuriousClones = getSpuriousSubclones(upper_bound, tp_cluster_prob)
    noOfiter = 1
    finalTree = Tree # initialize final tree
    noOfMutations = len(D_matrix[0])
    print("No of mutations ",noOfMutations)
    finalClusterCells = tp_cluster_cells
    finalClusterGen = tp_cluster_genotype
    print(" ==== Estimated FP FN before eliminating any spurious subclones ==== ")
    tp_FP, tp_FN = timepointFPFN(D_matrix, tp_cluster_cells, tp_cluster_genotype)
    print(" Estimated FP ",tp_FP," Estimated FN ",tp_FN)

    # Before saving the plots remove all exisitng plots including hidden files
    #if os.path.exists("plots/"+plotOp):
    #    items = os.listdir("plots/"+plotOp)
    #    for item in items:
    #        item_path = os.path.join("plots/"+plotOp, item)
    #        os.remove(item_path)
    #        if item.startswith('.'):
    #            os.remove(item_path)
    #else:
    #    os.makedirs("plots/"+plotOp) # Make the dir to save the plots

    #for tp_cluster in tp_cluster_prob:
    for tp_cluster in spuriousClones:
        print(" Iterations ============================== ",noOfiter)
        print("TP cluster gen ",tp_cluster_genotype)
        print("TP cluster cells ",tp_cluster_cells)
        print("Tp_cluster_prob ",tp_cluster_prob)
        # Update these for every iteration
        #new_tp_cluster_genotype = tp_cluster_genotype
        #new_tp_cluster_cells = tp_cluster_cells
        #print(tp_cluster)
        #if tp_cluster not in tp_cluster_prob:
        #    continue
        tp = int(tp_cluster.split('_')[0])
        cluster = int(tp_cluster.split('_')[1])
        print(" Timepoint ",tp," ignore cluster ",cluster)
        # save the cluster genotype before correction to use it later
        tp_cg_beforeCorr = copy.deepcopy(tp_cluster_genotype)
        tp_cc_beforeCorr = copy.deepcopy(tp_cluster_cells)

        cells = tp_cluster_cells[tp][cluster]
        alpha = tp_alpha[tp]
        beta = tp_beta[tp]
        #print("Cluster genotype ",tp_cluster_genotype[tp])
        # Before any changes in getting the Tree by removing the spurious subclones note the FP FN
        old_tp_FP, old_tp_FN = timepointFPFN(D_matrix, tp_cluster_cells, tp_cluster_genotype)
        print(" Iteration ",noOfiter," old FP ",old_tp_FP," old FN ",old_tp_FN)

        # Reassign cells and update the cluster genotype after removing spurious subclone
        cell_cluster, cluster_genotype = cell_cluster_likelihood(D_matrix, alpha, beta, cells, tp_cluster_genotype[tp], cluster, tp_cluster_cells[tp])
        # If nothing to reassign then we break with Tree with highest prob
        if cell_cluster == "Done" and cluster_genotype == "Done":
            break

        tp_cluster_cells[tp] = cell_cluster
        # Update new cluster genotypes
        tp_cluster_genotype[tp] = cluster_genotype
        print(" Cluster Genotype before reassigning cells =============== ")
        printMutFromGenotype(tp_cluster_genotype[tp])
        # After cell reassigning get the updated corrected cluster genotype
        tp_cluster_cells, tp_cluster_genotype = correctClustering(D_matrix, tp_cluster_cells, tp_alpha, tp_beta)
        print("Cluster Genotype after reassigning cells =============== ")
        printMutFromGenotype(tp_cluster_genotype[tp])

        # Prob after merging clones and reassigning cells
        afterMergeProb = calc_cluster_prob(tp_cluster_cells, tp_cluster_genotype, D_matrix, tp_alpha, tp_beta, tp_MR)

        # Get the Tree 
        Tree, clone_node = getTree(tp_cluster_cells, tp_cluster_genotype)
        print(" Before correction iteration ",noOfiter," Clone node ",clone_node)
        print("Tp_cluster_prob ",tp_cluster_prob)
        print(" After merging prob ",afterMergeProb)
        
        # Fix the parallel and back mutation. After correction the tp_cluster_genotype will change. 
        parallel_mut, parallelMut_edges = getParallelMut(Tree)
        backMut_edges = getBackMutCount(Tree)
        
        # Plot tree before parallel and back mutation correction
        plotTree(Tree, parallelMut_edges, backMut_edges,"plots/"+plotOp+"/iter"+str(noOfiter)+"_bc.png", "Iter "+str(noOfiter)+" before correction", sample)

        print(" Parallel edges ",parallelMut_edges)
        parallelEdges = finalParallelMutEdges(Tree, parallel_mut, parallelMut_edges, D_matrix, tp_alpha, tp_beta)
        backEdges = finalBackMutEdges(Tree, backMut_edges, D_matrix, tp_alpha, tp_beta, k)
        print("Parallel edges ",parallelEdges)
        print("Back edges ",backEdges)
        # Update this to keep only back mutations when there is an overlap of parallel and back mutations
        #parallelEdges = filterParallelEdges(parallelEdges, backEdges)

        Tree = correctParallelMut(Tree, parallelMut_edges, parallelEdges)
        printTree(Tree)
        Tree = correctBackMut(Tree, backMut_edges, backEdges)
        Tree = recheckTree(Tree)
        print("========== TREE after correcting PARALLEL and BACK MUTATIONS =========== ")
        printTree(Tree)

        # Genotype after correcting parallel and back mutations
        new_tp_cluster_genotype = genotypeFromTree(Tree, tp_cluster_genotype, clone_node, noOfMutations)

        #new_backMut_edges = getBackMutCount(Tree)
        #new_plMut, new_plMut_edges = getParallelMut(Tree)
        #plotTree(Tree, new_plMut_edges, new_backMut_edges,"plots/"+plotOp+"/iter"+str(noOfiter)+"_ac.png", "Iter "+str(noOfiter)+" after correction")

        # After updating the Tree get set of new FP FN
        new_tp_FP, new_tp_FN = timepointFPFN(D_matrix, tp_cluster_cells, new_tp_cluster_genotype)
        print(new_tp_cluster_genotype)
        print("=========== DIFF FP FN ===============")
        tp_diffFP, tp_diffFN = diffFPFN(old_tp_FP, old_tp_FN, new_tp_FP, new_tp_FN)
        #if userRate > 0:
        #    dropSubclone = checkDiff_FPFN(tp_diffFP, tp_diffFN, tp, userRate)
        #    if dropSubclone: # Reverse reassigning cells and prevent this subclone from dropping
        #        print(" Timepoint cluster cells before deciding to drop cluster ",tp_cluster_cells)
        #        tp_cluster_cells = tp_cc_beforeCorr
        #        print(" Timepoint cluster cells after deciding to drop cluster ",tp_cluster_cells)
        #        print(" Timepoint cluster genotype before deciding to drop cluster ",tp_cluster_genotype)
                # Update tp_cluster_genotype to have previous value
        #        tp_cluster_genotype = tp_cg_beforeCorr
        #        print(" Timepoint cluster genotype after not dropping the cluster ",tp_cluster_genotype)
        #        print("Iteration ",noOfiter," Cluster NOT DROPPED ========= ",cluster)

        #        noOfiter = noOfiter+1
        #        continue

        new_backMut_edges = getBackMutCount(Tree)
        new_plMut, new_plMut_edges = getParallelMut(Tree)
        plotTree(Tree, new_plMut_edges, new_backMut_edges,"plots/"+plotOp+"/iter"+str(noOfiter)+"_ac.png", "Iter "+str(noOfiter)+" after correction", sample)

        #print(" Cells in each TP cluster ",tp_cluster_cells)
        
        new_tp_cluster_prob = calc_cluster_prob(tp_cluster_cells, new_tp_cluster_genotype, D_matrix, tp_alpha, tp_beta, tp_MR)
        print("After correction Iteration ",noOfiter," Clone node ",clone_node)
        print("Tp_cluster_prob ",new_tp_cluster_prob)

        new_tree_prob = calculate_treeProb(new_tp_cluster_prob)
        print("Tree prob ",tree_prob," New tree prob ",new_tree_prob)
        noOfiter = noOfiter+1
        
        if new_tree_prob > tree_prob:
            tree_prob = new_tree_prob
            finalTree = Tree
            finalClusterCells = tp_cluster_cells
            finalClusterGen = new_tp_cluster_genotype
            final_backMut_edges = new_backMut_edges
            final_plMut_edges = new_plMut_edges
            print("FINAL TREE ITERATION ",noOfiter-1)
            #print("Back mutations ",final_backMut_edges)
            # plot the final tree. Last one gets overwritten
            plotTree(finalTree, final_plMut_edges, final_backMut_edges,"plots/"+plotOp+"/"+bnpcRun+"_Final_bc.png", "Final Tree",sample)
        else:
            print("=========== TREE NOT SELECTED =============")
            printTree(Tree)
            continue
    # the saved Tree will be the one with highest prob
    print(" ================= FINAL TREE =================== ")
    finalTree = recheckFinalTree(finalTree)
    printTree(finalTree)
    plotTree(finalTree, final_plMut_edges, final_backMut_edges,"plots/"+plotOp+"/"+bnpcRun+"_Final_ac.png", "Final Tree ", sample)
    print("Final tree prob ",tree_prob)
    return finalTree, final_backMut_edges, tree_prob

# Input is cells in each timepoint file.
# Output a dictionary with timepoint as key and cells as value.
def cellTimepoints(tpCellsFile):
    tcf = open(tpCellsFile, 'r')
    tcf_lines = tcf.readline().rstrip('\n')
    tpCells = []
    while(tcf_lines != ""):
        #tp = (tcf_lines.split('\t')[0]).split('_')[1]
        cells = (tcf_lines.split('\t')[1]).split(';')
        print(cells)
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

''' Get the missing rate for each time point. '''
# Input is cells in each time point and noisy D matrix
def calcFPFN(D,G):
    FN = 0
    FP = 0
    TP = 0
    TN = 0
    #print(" Mut ",len(D[0]))
    #print(" Mut ",len(G[0]))
    for i in range(len(G)):
        for j in range(len(G[i])):
            if D[i][j] == 0 and G[i][j] == 1:
                FN = FN + 1
            if D[i][j] == 1 and G[i][j] == 0:
                FP = FP + 1
            if D[i][j] == 1 and G[i][j] == 1:
                TP = TP + 1
            if D[i][j] == 0 and G[i][j] == 0:
                TN = TN + 1
    if FP+TN == 0:
        FPR = 0
    else:
        FPR = FP/(FP+TN)

    if FN+TP == 0:
        FNR = 0
    else:
        FNR = FN/(FN+TP)
    return FPR, FNR

''' Calculate estimated FP FN rates for each timepoint after getting the final tree. '''
def timepointFPFN(D_matrix, finalClusterCells, finalClusterGen):
    tp_FP = {}
    tp_FN = {}
    for tp, cluster in finalClusterCells.items():
        tp_cluster_G = []
        tp_D = []
        tp_cells = []

        for c in cluster:
            cells = finalClusterCells[tp][c]
            tp_cells.extend(cells)
            for cell in cells:
                tp_cluster_G.append(finalClusterGen[tp][c])
       
        print("Tp cluster cells ",len(tp_cells))
        #print(len(tp_cluster_G))
        for c in tp_cells:
            tp_D.append(D_matrix[c])
        #print(len(tp_D))

        FPR, FNR = calcFPFN(tp_D, tp_cluster_G)
        tp_FP[tp] = FPR
        tp_FN[tp] = FNR
        #print(" Timepoint ",tp," FP ",FPR," FN ",FNR)
    return tp_FP, tp_FN

''' Get the initial clustering results. '''
def intialClusterResults(tpClusters, allClusters, tpCells, D_matrix, tp_alpha, tp_beta, tp_MR):
    tp_cluster_cells = getNewClusters(tpClusters, allClusters, tpCells)
    tp_reassignedCells, tp_updatedCG = correctClustering(D_matrix, tp_cluster_cells, tp_alpha, tp_beta)
    tp_cluster_prob = calc_cluster_prob(tp_reassignedCells, tp_updatedCG, D_matrix, tp_alpha, tp_beta, tp_MR)
    print(" Timepoint cluster prob ",tp_cluster_prob)
    sorted_cluster_prob = dict(sorted(tp_cluster_prob.items(), key=lambda item:item[1]))
    print(sorted_cluster_prob)
    return sorted_cluster_prob, tp_reassignedCells, tp_updatedCG

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
#parser.add_argument("-k", "--k", help="No. of losses allowed")
#parser.add_argument("-e", "--e", help="FP FN error rate expected")
#parser.add_argument("-op", "--op", help="Path to save the resulting Tree.")
#parser.add_argument("-sample", "--sample", help="Sample name to plot the graphs.")
# User should be able to enter FP and FN rates if known while running BnpC and enter those values here.
#args = parser.parse_args()

#tpCells = cellTimepoints(args.cells)
#D_matrix = readDMatrix(args.D)
#tp_alpha, tp_beta = readFPFNvalues(args.logFile)
#tp_MR = timepoint_missingRate(tpCells, D_matrix)

# First iteration before eliminating spurious subclones 
# ======== Clustering results and calculating their probabilities ==========
#tp_cluster_cells = getNewClusters(args.tclusters, args.cluster, tpCells)
#tp_reassignedCells, tp_updatedCG = correctClustering(D_matrix, tp_cluster_cells, tp_alpha, tp_beta)

#tp_cluster_prob = calc_cluster_prob(tp_reassignedCells, tp_updatedCG, D_matrix, tp_alpha, tp_beta, tp_MR)
#print(" Timepoint cluster prob ",tp_cluster_prob)
#sorted_cluster_prob = dict(sorted(tp_cluster_prob.items(), key=lambda item:item[1]))
#print(sorted_cluster_prob)


#Tree, finalClusterCells, finalClusterGen = selectOptimalTree(D_matrix, tp_alpha, tp_beta, sorted_cluster_prob, tp_reassignedCells, tp_updatedCG, tpCells, int(args.k), float(args.e), args.op, args.sample)
#print(" ================= FINAL TREE =================== ")
#printTree(Tree)
#print(finalClusterCells)
#print(finalClusterGen)

#finalClusterCells[0] = {1: [6, 8, 9, 21, 36, 42, 43, 51, 56, 59, 63, 65, 57, 58, 23, 37, 74, 85, 5, 15, 71, 79, 86, 88, 12, 13, 19]}
#finalClusterGen[0] = {1: [1, 0, 1, 1, 1, 1, 1, 0, 0, 0, 1, 0, 0, 1, 0, 1, 1, 1, 1, 1, 0, 0, 1, 1, 0, 0, 1, 1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 0, 0, 1, 1, 1, 0, 1, 1, 0, 0, 0, 1, 1, 1, 1, 1]}
#timepointFPFN(D_matrix, finalClusterCells, finalClusterGen)

