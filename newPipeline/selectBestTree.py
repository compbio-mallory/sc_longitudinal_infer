import argparse
import os
import math
import numpy as np
from selectOptimalTree import cellTimepoints, readDMatrix, timepoint_missingRate, readFPFNvalues, intialClusterResults, selectOptimalTree, plotTree, saveTree

# Calculate all possible list of FP or FN rates for each timepoint given the FP and FN for each timepoint from the BnpC runs
# The list starts from median and ends with max value. It has all range of values including step size
# Input is an array of FP OR FN rates for each timepoint from all BnpC runs and step size for the errors 
def possibleListOfErrors(all_bnpc_errors, step_size):
    tp_error = {}
    for i in range(len(all_bnpc_errors)):
        maxVal = np.max(all_bnpc_errors[i])
        median = np.median(all_bnpc_errors[i])
        print("Median ",median," maxVal ",maxVal)
        errors = []
        errors.append(median)
        value = median
        while(value <= maxVal):
            print("Step size ",step_size)
            value = value + step_size
            print("Value ",value)
            errors.append(value)
            #step_size = step_size + step_size
        if median != maxVal:
            errors.append(maxVal)
        tp_error[i] = errors
    return tp_error

# Get the FP and FN rates for each timepoint from all BnpC runs and gets the thresholds based on the step size provided.
# Input is no. of bnpc runs, no. of timepoints, clustering results location
def getFPFN_threshold(m, t, cLoc):
    all_bnpc_FP = []
    all_bnpc_FN = []
    epsilon = float(1e-6)
    # create empty lists for each timepoint to save the FP FN values.
    for j in range(len(t)):
        all_bnpc_FP.append([])
        all_bnpc_FN.append([])
    
    for i in range(1,m+1):
        for j in range(len(t)):
            # Read the alpha and beta values here for each timepoint.
            # Then assign them to all_bnpc_FP[j], all_bnpc_FN[j]
            lfile = open(cLoc+"/m"+str(i)+"/"+t[j]+"/errors.txt",'r')
            lf_line = lfile.readline().rstrip('\n')
            lf_line = lfile.readline().rstrip('\n')
            lf_arr = lf_line.split('\t')
            fn_rate = float(lf_arr[3]) + epsilon
            fp_rate = float(lf_arr[5]) + epsilon
            all_bnpc_FP[j].append(fp_rate)
            all_bnpc_FN[j].append(fn_rate)
    print("All bnpc FP ",all_bnpc_FP)
    print("All bnpc FN ",all_bnpc_FN)
    tp_FP = possibleListOfErrors(all_bnpc_FP, 0.001)
    tp_FN = possibleListOfErrors(all_bnpc_FN, 0.02)
    print(tp_FP," ",tp_FN)
    return tp_FP, tp_FN

#Function to create a recursive for loop because we will have different no. of timepoints and varying no. of error rates
# Output is list having inner lists having pair of FP and FN rates for each timepoint. Ex: 1 combination for 3 timepoints: [[[FP0 , FN0], [FP1, FN1], [FP2, FN2]]]
def for_recursive(number_of_loops, resultList, range_list, current_index=0, iter_list = []):
  if iter_list == []:
    iter_list = [0]*number_of_loops

  if current_index == number_of_loops-1:
    for iter_list[current_index] in range_list[current_index]:
        resultList.append(iter_list.copy())
  else:
    for iter_list[current_index] in range_list[current_index]:
        for_recursive(number_of_loops, resultList = resultList, iter_list = iter_list, range_list = range_list,  current_index = current_index+1) 
  return resultList

#Create a list of list of dictionaries having possible combinations of FP FN rates per timepoint.
# Ex. [{0: 0.001, 1: 0.01, 2: 0.05}, {0: 0.1, 1: 0.2, 2: 0.3}]
def combinationsOfFPFNs(tp_FP, tp_FN):
    FP_list = []
    FN_list = []
    tp_pair_errors = {}
    # Make separate lists for FP,FN for each timepoint.
    tp_keys = list(tp_FP.keys())
    tp_values = list(tp_FP.values())

    for tp in tp_FP:
        FPrates = tp_FP[tp]
        FNrates = tp_FN[tp]
        for fp in FPrates:
            for fn in FNrates:
                pair_errors = []
                pair_errors.append(fp) # 0 index has the FPs
                pair_errors.append(fn) # 1 index has the FNs
                if tp in tp_pair_errors:
                    tp_pair_errors[tp].append(pair_errors)
                else:
                    tp_pair_errors[tp] = [pair_errors]
    # Sample tp_pair_errors. key -> timepoint, values: [possible combinations of FP and FN rates]
    #tp_pair_errors = {0: [[0.007993225, 0.11400099999999999], [0.007993225, 0.13400099999999998], [0.007993225, 0.137301]], 1: [[0.02069066, 0.093001],[0.03, 0.09]], 2: [[1e-06, 0.064301],[0.04, 0.1]]}
    print("All possible timepoint FP FNs ",tp_pair_errors)
    # Using recursion to get all possible combination of FP FNs for each timepoint. Since no. of timepoints may vary we use dynamic no. of for loops
    # Ex: 1 combinationsOfErrors for 3 timepoints: [[[FP0 , FN0], [FP1, FN1], [FP2, FN2]]]
    combinationsOfErrors = for_recursive(range_list = list(tp_pair_errors.values()), resultList = [], number_of_loops=len(tp_keys))
    print("Return from recursion ",combinationsOfErrors)
    return combinationsOfErrors

''' Return the Tree with the highest prob from all the BnpC runs. '''
# Input: m = no. of Bnpc runs, t = no. of timepoints, cLoc = clustering results location, cells in each timepoint, D matrix, no. of losses allowed, FP FN rates, path to save plots, sample to visualize the plot
def getTreeWithHighestProb(m, t, cLoc, tp_MR, tpCells, D_matrix, k, combinationsOfErrors, plotOp, sample):
    Tree_prob = {}
    Tree_backMut = {}
    # Save the estimated error rates by the tree with the max prob
    bnpc_FP = {}
    bnpc_FN = {}
    # Save the thresholds used to get the Tree with max prob
    FP_threshold = {}
    FN_threshold = {}
    for i in range(1,m+1):
    #for i in range(1,2):
        #allClusters = cLoc+"/m"+str(i)+"/all/assignment.txt" # assignment.txt from BnpC run across all timepoints
        print("BnpC run ",i)
        #Tree_prob = {}
        #FP_parallel_back = {}
        #FN_parallel_back = {}
        #FP_threshold = {}
        #FN_threshold = {}

        tp_files = [] # save timepoint assignment.txt files
        tp_errors = [] # save timepoint FP FN rates
        for j in t:
            tp_files.append(cLoc+"/m"+str(i)+"/"+j+"/assignment.txt")
            tp_errors.append(cLoc+"/m"+str(i)+"/"+j+"/errors.txt")
        print("Timepoint file ",tp_files)
        tp_alpha, tp_beta = readFPFNvalues(tp_errors)
        #fp = max(tp_alpha.values()) / 0.8
        #fn = max(tp_beta.values()) / 0.8
        #print("Max TP alpha ",fp)
        #print("Max TP beta ",fn)
        #sorted_cluster_prob, tp_reassignedCells, tp_updatedCG = intialClusterResults(tp_files, allClusters, tpCells, D_matrix, tp_alpha, tp_beta, tp_MR)
        #sorted_cluster_prob, tp_reassignedCells, tp_updatedCG = intialClusterResults(tp_files, "", tpCells, D_matrix, tp_alpha, tp_beta, tp_MR)

        # Grid search with all possible combinations of FP and FN rates from combinationsOfErrors
        for ce in range(len(combinationsOfErrors)):
            tp_fp = {}
            tp_fn = {}
            for j in range(len(combinationsOfErrors[ce])):
                tp_fp[j] = combinationsOfErrors[ce][j][0]
                tp_fn[j] = combinationsOfErrors[ce][j][1]
            print("Grid search FPs ",tp_fp)
            print("Grid search FNs ",tp_fn)
            try:
                sorted_cluster_prob, tp_reassignedCells, tp_updatedCG = intialClusterResults(tp_files, "", tpCells, D_matrix, tp_alpha, tp_beta, tp_MR)
                Tree, back_mut, tree_prob = selectOptimalTree(D_matrix, tp_alpha, tp_beta, tp_MR, sorted_cluster_prob, tp_reassignedCells, tp_updatedCG, tpCells, k, tp_fp, tp_fn, plotOp, str(i), sample)
                bnpc_FP[tree_prob.copy()] = tp_alpha
                bnpc_FN[tree_prob.copy()] = tp_beta
                FP_threshold[tree_prob.copy()] = tp_fp
                FN_threshold[tree_prob.copy()] = tp_fn
                Tree_prob[tree_prob.copy()] = Tree.copy()
                Tree_backMut[tree_prob.copy()] = back_mut.copy()
                #FP_parallel_back[tree_prob] = final_tp_FP
                #FN_parallel_back[tree_prob] = final_tp_FN
            except:
                # some combination of errors might throw error while reassigning cells to the clusters because of mismatched error rates and infer wrong Tree hence we reject those Trees and continue with next combinations
                continue

        # Can use the following lines to save results per BnpC run
        #max_prob = max(Tree_prob)
        #print(list(Tree_prob.keys()))
        #print("Max prob ",max_prob)
        #finalTree = Tree_prob[max_prob]
        #finalBackMut = Tree_backMut[max_prob]
        #print("BnpC run ",i," BnpC FP rates ",tp_alpha," BnpC FN rates ",tp_beta)
        #print("FP thresholds used ",FP_threshold[max_prob]," FN thresholds used ",FN_threshold[max_prob])
        #print("FP FN rate used to fix parallel and back mutations ",FP_parallel_back," ",FN_parallel_back)
        #plotTree will save the .pdf file ans saveTree will save the Tree in .csv 
        #plotTree(finalTree, {}, finalBackMut,"plots/"+plotOp+"_"+str(i)+".pdf", "", sample)

    max_prob = max(Tree_prob)
    print(list(Tree_prob.keys()))
    print("Max prob ",max_prob)
    finalTree = Tree_prob[max_prob]
    finalBackMut = Tree_backMut[max_prob]
    print("Max FP rates ",bnpc_FP[max_prob]," Max FN rates ",bnpc_FN[max_prob])
    print("FP thresholds used ",FP_threshold[max_prob]," FN thresholds used ",FN_threshold[max_prob])
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
#parser.add_argument("-FP", "--FP", help="FP error rate expected", default=0.05)
#parser.add_argument("-FN", "--FN", help="FN error rate expected", default=0.35)
parser.add_argument("-op", "--op", help="Path to save the resulting Tree.")
parser.add_argument("-sample", "--sample", help="Sample name to plot the graphs.", default="default")
args = parser.parse_args()

tpCells = cellTimepoints(args.cells)
D_matrix = readDMatrix(args.D)
tp_MR = timepoint_missingRate(tpCells, D_matrix)
tp_FP, tp_FN = getFPFN_threshold(int(args.m), args.t, args.loc)
combinationsOfErrors = combinationsOfFPFNs(tp_FP, tp_FN)
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

getTreeWithHighestProb(int(args.m), args.t, args.loc, tp_MR, tpCells, D_matrix, int(args.k), combinationsOfErrors, args.op, args.sample)

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
