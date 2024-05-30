import argparse
import os
from selectOptimalTree import cellTimepoints, readDMatrix, timepoint_missingRate, readFPFNvalues, intialClusterResults, selectOptimalTree, plotTree, saveTree

''' Return the Tree with the highest prob from all the BnpC runs. '''
# Input: m = no. of Bnpc runs, t = no. of timepoints, cLoc = clustering results location
def getTreeWithHighestProb(m, t, cLoc, tp_MR, tpCells, D_matrix, k, fp, fn, plotOp, sample):
    Tree_prob = {}
    Tree_backMut = {}
    # Save the estimated error rates by the tree with the highest prob
    bnpc_FP = {}
    bnpc_FN = {}
    for i in range(1,m+1):
    #for i in range(3,4):
        #allClusters = cLoc+"/m"+str(i)+"/all/assignment.txt" # assignment.txt from BnpC run across all timepoints
        print("BnpC run ",i)
        tp_files = [] # save timepoint assignment.txt files
        tp_errors = [] # save timepoint FP FN rates
        for j in t:
            tp_files.append(cLoc+"/m"+str(i)+"/"+j+"/assignment.txt")
            tp_errors.append(cLoc+"/m"+str(i)+"/"+j+"/errors.txt")
        print("Timepoint file ",tp_files)
        tp_alpha, tp_beta = readFPFNvalues(tp_errors)
        fp = max(tp_alpha.values()) / 0.8
        fn = max(tp_beta.values()) / 0.8
        print("Max TP alpha ",fp)
        print("Max TP beta ",fn)
        #sorted_cluster_prob, tp_reassignedCells, tp_updatedCG = intialClusterResults(tp_files, allClusters, tpCells, D_matrix, tp_alpha, tp_beta, tp_MR)
        sorted_cluster_prob, tp_reassignedCells, tp_updatedCG = intialClusterResults(tp_files, "", tpCells, D_matrix, tp_alpha, tp_beta, tp_MR)
        Tree, back_mut, tree_prob = selectOptimalTree(D_matrix, tp_alpha, tp_beta, tp_MR, sorted_cluster_prob, tp_reassignedCells, tp_updatedCG, tpCells, k, fp, fn, plotOp, str(i), sample)
        bnpc_FP[tree_prob] = tp_alpha
        bnpc_FN[tree_prob] = tp_beta
        Tree_prob[tree_prob] = Tree
        Tree_backMut[tree_prob] = back_mut
    max_prob = max(Tree_prob)
    print(list(Tree_prob.keys()))
    print("Max prob ",max_prob)
    finalTree = Tree_prob[max_prob]
    finalBackMut = Tree_backMut[max_prob]
    print("Max FP rates ",bnpc_FP[max_prob]," Max FN rates ",bnpc_FN[max_prob])
    #plotTree(finalTree, {}, finalBackMut,"plots/"+plotOp+"/smallSA501.pdf", "", sample)
    # plotTree will save the .pdf file ans saveTree will save the Tree in .csv 
    plotTree(finalTree, {}, finalBackMut,"plots/"+plotOp+".pdf", "", sample)
    saveTree(finalTree, "TreeCSVs/"+plotOp+".csv")

parser = argparse.ArgumentParser()
#parser.add_argument("-tp", "--tp",dest = "tp", help="Timepoint cells genotype.")
#parser.add_argument("-cg", "--cg",dest = "cg", help="Consensus genotype from BnpC.")
# Get the timepoint clusters input as an array because we won't know how many timepoints we will have.
# Get the cells assigned in each timepoint instead of using index.
parser.add_argument("-m", "--m", help="No. of BnpC runs.")
parser.add_argument("-t", "--t", nargs="*", help="No. of time points.")
parser.add_argument("-loc", "--loc", help="Location to save the clustering results.")
parser.add_argument("-cells","--cells", help="Cells assigned at each timepoint.")
#parser.add_argument("-cluster", "--cluster",dest="cluster", help="Assignment.txt or file having cluster information")
#parser.add_argument("-f", "--f",dest="f", help="File having FP FN rates estimated from clusters")
parser.add_argument("-D","--D", help="D matrix.")
parser.add_argument("-k", "--k", help="No. of losses allowed")
parser.add_argument("-FP", "--FP", help="FP error rate expected", default=0.05)
parser.add_argument("-FN", "--FN", help="FN error rate expected", default=0.35)
parser.add_argument("-op", "--op", help="Path to save the resulting Tree.")
parser.add_argument("-sample", "--sample", help="Sample name to plot the graphs.", default="default")
args = parser.parse_args()

tpCells = cellTimepoints(args.cells)
D_matrix = readDMatrix(args.D)
tp_MR = timepoint_missingRate(tpCells, D_matrix)

# Clear the plots saving directory before saving the plots
#if os.path.exists("plots/"+args.op):
#    items = os.listdir("plots/"+args.op)
#    for item in items:
#        item_path = os.path.join("plots/"+args.op, item)
#        os.remove(item_path)
#        if item.startswith('.'):
#            os.remove(item_path)
#else:
#    os.makedirs("plots/"+args.op) # Make the dir to save the plots

getTreeWithHighestProb(int(args.m), args.t, args.loc, tp_MR, tpCells, D_matrix, int(args.k), float(args.FP), float(args.FN), args.op, args.sample)

# -m 10 -t X1 X2 X4 -loc largeSA501_bnpc -cells SA501/largerData_1/SA501.cell_timepoints.csv -D SA501/largerData_1/SA501.input.D.csv -k 0 -e 0 -op largeSA501_multiBnpC -sample large 
# Get the time point clusters folder and no. of runs as an input
# For each BnpC run read their assignment.txt from across all time point runs and their respective timepoint runs
#noOfRuns = args.m
#noOfTimepoints = args.t

#Tree = {} # key will be prob and value will be the saved Tree
#for i in range(m):
    # Get the BnpC across all time points assignment.txt file
#    for j in range(n):
        # Get the time point runs assignment.txt file
#        new_Tree, new_tree_prob = selectOptimalTree
#        if new_tree_prob > tree_prob:
#            tree_prob = new_tree_prob
#            Tree = new_Tree

# Resulting Tree is the tree with the highest prob across all runs
