"""
From this file first we check the parallel mutations for each runs.
We select the clustering results with min parallel mutations.
If there are multiple clustering results that have same least number of
parallel mutations then select one with least number of clusters.
Otherwise choose one randomly.
"""
import argparse
import os
import subprocess

''' Reads the bnpc assignment.txt file to calculate the no of clusters and returns the run with min cluster no. '''
def choose_best_cluster(min_mut_cluster, bnpc_path): # Read from BnpC assignment file
   #print(" BnpC Path ",bnpc_path)
   runNo_clusters = {}
   for run_no in min_mut_cluster:
       assignment_file = open(bnpc_path+"/"+run_no+"/assignment.txt",'r')
       af_line = assignment_file.readline().rstrip("\n")
       l_count = 0
       while(af_line != ""):
           if l_count == 0:
               l_count = l_count+1
               af_line = assignment_file.readline().rstrip("\n")
               continue
           clusters = (af_line.split('\t')[2]).split(' ')
           runNo_clusters[run_no] = len(set(clusters))
           af_line = assignment_file.readline().rstrip("\n")
   print(" Cluster nos for each run ",runNo_clusters)
   min_runNo = min(runNo_clusters, key=runNo_clusters.get)
   print(" Run ",min_runNo)
   return min_runNo

''' Finds the tree or trees with minimum number of parallel mutations. '''
def calculateParallelMutations(treeFiles, bnpc_path):
    # While looping through the file make a mutation list with the mutations added.
    # make a dictionary with node as key and mutations as list in value.
    # Loop over the mutation list and dict to check if mutation appears. 
    # While doing the above step maintain a dictionary where mutation is key and count as value.
    # Now find the min count of this dictionary value and save that in a dictionary with run number as key and min count value.
    # Then, finally select the cluster file with min count value.
    # If more than 1 cluster file with same min count value then check for no of clusters.
    cluster_parallel_mut = {} # Save the no. of parallel mutations
    for treef in treeFiles:
        cluster_name = (treef.split('/')[1]).split('_')[2]
        #print(" Cluster name ",cluster_name)
        tf = open(treef,'r')
        tf_line = tf.readline().rstrip("\n")
        mut_list = []
        node_mut = {}
        child_parent = {}
        while(tf_line != ""):
            node_id = tf_line.split('\t')[0]
            if node_id == '0': # Skipping the root node here
                tf_line = tf.readline().rstrip("\n")
                continue
            timepoint = tf_line.split('\t')[1]
            pid = tf_line.split('\t')[2]
            mut = (tf_line.split('\t')[4]).split(',')
            mut_list.extend(mut)
            node_mut[node_id] = mut
            child_parent[node_id] = pid
            tf_line = tf.readline().rstrip("\n")
        print(node_mut)
        mut_list = set(mut_list)
        #print(" Mutations ",mut_list)
        
        updated_node_mut = {}
        for child, parent in child_parent.items(): # Find the actual mutation set at a timepoint for each node.
            if parent == '0':
                updated_node_mut[child] = node_mut[child]
                continue
            if "," in parent:
                split_parent = parent.split(",")
                p_mut = node_mut[split_parent[0]]
            else:
                p_mut = node_mut[parent]
            child_mut = node_mut[child]
            actual_mutation = set(child_mut) - set(p_mut) # Actual mutation at a timepoint
            print(" Actual mutation of ",child," ",list(actual_mutation))
            updated_node_mut[child] = list(actual_mutation)

        print(" Updated node mutation ",updated_node_mut)
        mut_count_dict = {}
        for mutation in mut_list:
            mut_count = 0
            for n, m in updated_node_mut.items():
                if mutation in m:
                    mut_count = mut_count+1
            mut_count_dict[mutation] = mut_count
        #print(" Mutation count ",mut_count_dict)

        parallel_mut_count = 0
        for m,c in mut_count_dict.items(): # Find the parallel mutation values for each tree.
            if c >= 2:
                parallel_mut_count = parallel_mut_count+1
        cluster_parallel_mut[cluster_name] = parallel_mut_count
        
    print(" Parallel mutation count ",cluster_parallel_mut)
    min_mut_count_key = min(cluster_parallel_mut, key=cluster_parallel_mut.get) # First find the key with min parallel mutations.
    min_mut_count = cluster_parallel_mut[min_mut_count_key] # Then the least parallel mutation value.
    min_mut_cluster = []
    for c,m in cluster_parallel_mut.items():
        if m == min_mut_count:
            min_mut_cluster.append(c)
    print(" Clusters with least parallel mutations ",min_mut_cluster)

    if len(min_mut_cluster) > 1: # Select the one with min no of clusters if more than one runs exist with least parallel mutations
        best_cluster_run = choose_best_cluster(min_mut_cluster, bnpc_path)
        base_path = (treeFiles[0].split('/')[1]).split('_')
        print(" Base path ",base_path)
        best_tree_path = base_path[0]+'_'+base_path[1]+'_'+best_cluster_run+'_tree.csv'
        print(" Best tree path ",best_tree_path)
        cp_cmd = 'cp '+treeFiles[0].split('/')[0]+"/"+best_tree_path+' bestTree/'
        subprocess.call(cp_cmd, shell=True)
        cp_cloneNodes = 'cp all_trees/'+base_path[0]+'_'+base_path[1]+'_'+best_cluster_run+'_cloneNode.npy'+' bestTree/'
        cp_cloneCells = 'cp all_trees/'+base_path[0]+'_'+base_path[1]+'_'+best_cluster_run+'_cloneCells.npy'+' bestTree/'
        subprocess.call(cp_cloneNodes, shell=True)
        subprocess.call(cp_cloneCells, shell=True)
    else:
        base_path = (treeFiles[0].split('/')[1]).split('_')
        best_tree_path = base_path[0]+'_'+base_path[1]+'_'+min_mut_cluster[0]+'_tree.csv'
        print(" Best tree path ",best_tree_path)
        cp_cmd = 'cp '+treeFiles[0].split('/')[0]+"/"+best_tree_path+' bestTree/'
        subprocess.call(cp_cmd, shell=True)
        cp_cloneNodes = 'cp all_trees/'+base_path[0]+'_'+base_path[1]+'_'+min_mut_cluster[0]+'_cloneNode.npy'+' bestTree/'
        cp_cloneCells = 'cp all_trees/'+base_path[0]+'_'+base_path[1]+'_'+min_mut_cluster[0]+'_cloneCells.npy'+' bestTree/'
        subprocess.call(cp_cloneNodes, shell=True)
        subprocess.call(cp_cloneCells, shell=True)

parser = argparse.ArgumentParser()
parser.add_argument("-tree", "--tree",dest ="tree", nargs="*", help="Final Longitudinal trees")
parser.add_argument("-bnpcOp", "--bnpcOp",dest="bnpcOp", help="BnpC output dir")
args = parser.parse_args()

print(" List of tree files ",args.tree)
calculateParallelMutations(args.tree, args.bnpcOp)
