import ete3
from ete3 import Tree
import argparse
import pandas as pd
import ast

''' Read the clonal phylogeny and return a dictionary of cells belonging to the clusters. '''
def readClonalPhylogeny(clonal_phylogeny):
    #t = Tree('less24hrs_t4/default/rep1/rep1/samples/best/best_MAP_tree.txt', format=1)
    t = Tree(clonal_phylogeny, format=1)
    cluster_cells_dict = {}
    # iterate to get the children of each node
    for node in t.traverse("postorder"):
        # get the children of the nodes
        children = node.get_children()
        for child in children:
        #    print(node.name," -- " ,child.name)
            if node.name in cluster_cells_dict:
                temp = cluster_cells_dict[node.name]
                temp.append(child.name)
                cluster_cells_dict[node.name] = temp
            else:
                temp = []
                temp.append(child.name)
                cluster_cells_dict[node.name] = temp
    #print(cluster_cells_dict)
    return cluster_cells_dict

''' Read the clonal phylogeny and include the internal nodes to connect all nodes. '''
def treeWithInternalNodes(clonal_phylogeny):
    f = open(clonal_phylogeny, 'r')
    line = f.readline()
    line = line.rstrip()
    f.close()
    newick_str = '('+ line[:-1] + ')R:0;'
    #print(newick_str)
    tree = Tree(newick_str, format=1)
    #newickWithnames = add_internal_node_names(tree)
    #print(newickWithnames)
    tree = Tree(newick_str, format=1)
    par2child = {}
    child2par = {}
    #child2par[] = -1
    node_count = 0
    for node in tree.traverse('preorder'):
        if not node.is_leaf() and not node.name:
            node.name = 'IN' + str(node_count)
            node_count += 1
        if node.name not in par2child:
            if node.name.startswith('sc'):
                continue
            par2child[node.name] = []
            for child in node.children:
                if child.name.startswith('sc'):
                    continue
                if not child.name:
                    child.name = 'IN' +str(node_count)
                    node_count += 1
                par2child[node.name].append(child.name)
                child2par[child.name] = node.name
    return par2child, child2par

''' Get the clonal genotype by majority voting of the cells. '''
def getClonalGenotype(cells_genotype):
    clonal_genotype = []
    for j in range(0,len(cells_genotype[0])):
        total1s = 0
        total0s = 0
        for i in range(0,len(cells_genotype)):
            if cells_genotype[i][j] == 1:
                total1s = total1s+1
            else:
                total0s = total0s+1
        if total1s >= total0s:
            clonal_genotype.append(1)
        if total0s > total1s:
            clonal_genotype.append(0)
    return clonal_genotype

''' Get the mutation of the clones from the clonal genotype. '''
def getClonalMutation(cluster_cells_dict, cell_genotype):
    cluster_mutation = {}
    for cluster, cells in cluster_cells_dict.items():
        #print(cluster," ",cells)
        if cluster == "":
            continue
        cells_genotype = []
        for cell in cells:
            cell_no = cell.split('c')[1]
            cells_genotype.append(cell_genotype[int(cell_no)])
            #print(cell_no)
        #print(cells_genotype)
        #print("==========================")
        clonal_genotype = getClonalGenotype(cells_genotype)
        #print(clonal_genotype," ",len(clonal_genotype)," ",len(cell_genotype[0]))
        
        mut_list = []
        # Get the mutations from clonal_genotype
        for i in range(len(clonal_genotype)):
            if clonal_genotype[i] == 1:
                mut_list.append(i)
        cluster_mutation[cluster] = mut_list
    #print(cluster_mutation)
    return cluster_mutation

''' Get the mutations for all nodes including the internal nodes. '''
def getNodeMutation(par2child, child2par, cluster_mutation):
    node_mut = {}
    cluster_inNodes = {} # Internal nodes connected to the clones. We can then assign cluster mutations to these nodes.
    for child, par in child2par.items(): # Filter to get only the internal nodes associated with child
        internalNode_list = []
        if 'C' in child:
            internalNode_list.append(par)
            #par_nodes = par2child[par]
            #for pn in par_nodes:
            #    if 'IN' in pn:
            #        internalNode_list.append(pn)

            cluster_inNodes[child] = internalNode_list
    #print(" Cluster internal nodes ",cluster_inNodes)           
    
    # Create mutation dictation with internal nodes having same mutations as clusters. 
    # Next, create a dictionary mentioning internal nodes having more than one cluster and hence common mutations from all the clusters.

    for par, child in par2child.items():
        if par == 'R' or 'C' in par: # No need to check for root node because it will not have any mutations
            continue

        c_mut = []
        for c in child: # Loop over the connected children and retrieve their mutations
            if 'C' in c: # Get the mutations from the clones/clusters
                if c_mut == []:
                    c_mut.extend(cluster_mutation[c])
                else:
                    temp_c_mut = cluster_mutation[c]
                    c_mut = list(set(temp_c_mut) & set(c_mut))

        cluster_mutation[par] = c_mut
    #print(" Cluster mutations ",cluster_mutation)
    
    # Once we have mutations for the internal nodes from the clusters we need to filter them based on their parents mutations.
    for par, child in par2child.items():
        if 'C' in par or par == 'R': # No need to check for root node because it will not have any mutations
            continue
        #updated_mut = []
        # If two internal nodes are children.
        checkInNode1 = child[0]
        checkInNode2 = child[1]
        if 'IN' in checkInNode1 and 'IN' in checkInNode2:
            cluster_mutation[par] = list(set(cluster_mutation[child[0]]).union(set(cluster_mutation[child[1]])))
            continue

        for c in child:
            if 'C' in c:
                temp_mut = set(cluster_mutation[c]) - set(cluster_mutation[par]) # Only include new mutations
                if len(temp_mut) > 0:
                    cluster_mutation[c] = list(temp_mut)
            else:
                cluster_mutation[c] = cluster_mutation[par] # Internal nodes will have mutations of its parent internal node.

            #if updated_mut == []:
            #    updated_mut.extend(cluster_mutation[c])
            #else:
                #temp_c_mut = cluster_mutation[c]
                #updated_mut = list(set(temp_c_mut) & set(updated_mut))

            #cluster_mutation[c] = cluster_mutation[par] # Internal nodes will have mutations of its parent internal node.
    print(" Updated Cluster mutations ",cluster_mutation)
    return cluster_mutation

# Also, get mut_edge dict where mutation is the key and parent_child node as values.
# Finding, SNVs on same branches will be same here.
def getEdgeMutation(par2child,updated_cluster_mutation):
    edge_mut = {}
    parallel_edges = []
    for par, child in par2child.items():
        if 'C' in par:
            continue
        if 'R' in par:
            edge_mut_key = par+'_'+child[0]
            edge_mut_val = updated_cluster_mutation[child[0]]
            edge_mut[edge_mut_key] = edge_mut_val
            continue
        node_parallel_edges = []
        for c in child:
            edge_mut_key = par+'_'+c
            print(" Parent ",par," mutation ",updated_cluster_mutation[par])
            print(" Child ",c," mutation ",updated_cluster_mutation[c])
            edge_mut_val = set(updated_cluster_mutation[c]) - set(updated_cluster_mutation[par])
            #edge_mut_val = updated_cluster_mutation[c]
            print(edge_mut_val)
            edge_mut[edge_mut_key] = edge_mut_val
            node_parallel_edges.append(edge_mut_key)
        parallel_edges.append(node_parallel_edges)
    print(" Edge mutation ",edge_mut)
    print(" Parallel edges ",parallel_edges)

    mut_edge = {} # Maintain a dict with mutation as the key and edges as values.
    for edge,mut in edge_mut.items():
        for m in mut:
            if m in mut_edge:
                t_mut_edge = mut_edge[m]
                t_mut_edge.append(edge)
                mut_edge[m] = t_mut_edge
            else:
                t_mut_edge = []
                t_mut_edge.append(edge)
                mut_edge[m] = t_mut_edge
    #print(mut_edge)
    return mut_edge,edge_mut,parallel_edges

# if pair of SNVs already checked then prevent checking
def checkPairedSNVs(sameSNVs,mut1,mut2):
    for mut in sameSNVs:
        mut_list = mut.split('_')
        #print(mut1," ",mut2," ",mut_list)
        if str(mut1) in mut_list and str(mut2) in mut_list:
            return True

# Get the pair of SNVs that are on the same branch
def sameEdgeMutation(mut_edge):
    sameSNVs = {} # pair of SNVs appearing on same edges together
    for mut1, edge1 in mut_edge.items():
        for mut2, edge2 in mut_edge.items():
            edge_list = []
            #print(" Edge 1 ",edge1)
            #print(" Edge 2 ",edge2)
            if checkPairedSNVs(sameSNVs,mut1,mut2) or mut1 == mut2:
                continue
            for edge in edge1:
                if edge in edge2:
                    edge_list.append(edge)
            #print(" Edge list ",edge_list)
            if edge_list != []:
                sameSNVs[str(mut1)+'_'+str(mut2)] = set(edge_list)
                
    #print(sameSNVs.keys())
    return sameSNVs

# SNVs appearing in parallel egdes
def parallelSNVs(edge_mut, parallel_edges):
    inComSNVs = {} # incomparable pair of SNVs appearing in parallel edges
    #mutations = mut_edge.keys()
    for pe in parallel_edges:
        edge1 = pe[0]
        edge2 = pe[1]
        edge1_mut = edge_mut[edge1]
        edge2_mut = edge_mut[edge2]
        for mut1 in edge1_mut:
            for mut2 in edge2_mut:
                if checkPairedSNVs(inComSNVs,mut1,mut2) or mut1 == mut2:
                    continue
                inComSNVs[str(mut1)+'_'+str(mut2)] = pe
    #print(inComSNVs.keys())
    return inComSNVs

# Get the ancestral relation of SNV pairs.
def getAncestralRelation(par2child, updated_cluster_mutation):
    #inComSNV_pairs = list(inComSNVs.keys())
    #print(inComSNV_pairs)
    ancestralSNVs = {}
    for par, child in par2child.items():
        if par == 'R':
            continue
        par_mut = updated_cluster_mutation[par]
        for c in child:
            c_mut = updated_cluster_mutation[c]
            #print(par+'_'+c," Child mut ",c_mut)
            for m1 in par_mut:
                for m2 in c_mut:
                    if checkPairedSNVs(ancestralSNVs,m1,m2) or m1 == m2:
                        continue
                    snv_pair = str(m1)+'_'+str(m2)
                    #if snv_pair in inComSNV_pairs: # not include incomparable snvs
                    #    continue
                    if snv_pair in ancestralSNVs:
                        t_edge = ancestralSNVs[snv_pair]
                        t_edge.append(par+'_'+c)
                        ancestralSNVs[snv_pair] = t_edge
                    else:
                        t_edge = []
                        t_edge.append(par+'_'+c)
                        ancestralSNVs[snv_pair] = t_edge
    #print(" Ancestral SNVs ",ancestralSNVs.keys())
    return ancestralSNVs

''' Get the incomparable SNVs. '''
def getInComparableSNVs(updated_cluster_mutation, sameSNVs, ancestralSNVs):
    inCompSNVs = []
    sameSNVs_list = set(sameSNVs.keys())
    ancestralSNVs_list = set(ancestralSNVs.keys())
    print(" sameSNVs no. ",len(sameSNVs_list)," ancestralSNVs no. ",len(ancestralSNVs_list))
    for node1, mut1 in updated_cluster_mutation.items():
        for node2, mut2 in updated_cluster_mutation.items():
            for m1 in mut1:
                for m2 in mut2:
                    if checkPairedSNVs(inCompSNVs,str(m1),str(m2)) or str(m1) == str(m2) or (str(m1)+'_'+str(m2) in sameSNVs_list) or (str(m1)+'_'+str(m2) in ancestralSNVs_list) or (str(m2)+'_'+str(m1) in sameSNVs_list) or (str(m2)+'_'+str(m1) in ancestralSNVs_list):
                        continue
                    inCompSNVs.append(str(m1)+'_'+str(m2))
    print(" Incomparable SNVs ",len(set(inCompSNVs)))
    return set(inCompSNVs)

''' save the SNV pairs to a file. '''
def save_files(SNVset,fileName):
    with open(fileName,'w') as f:
        f.write(str(SNVset))

    # check if the saved file can be read
    #with open(fileName,'r') as f:
    #    my_set = ast.literal_eval(f.read())
    #print(my_set)

# iterate to get the parent of each node
#for node in t.traverse("levelorder"):
#    if not node.is_root():
#        print(node.name," parent --> ",node.up.name)

parser = argparse.ArgumentParser()
parser.add_argument("-tree", "--tree",dest ="tree", help="Tree generated by SiCloneFit")
parser.add_argument("-genotype", "--genotype",dest ="genotype", help="Genotype saved by SiCloneFit")
parser.add_argument("-output", "--output",dest ="output", help="Directory to save the SNV pairs")
args = parser.parse_args()
cluster_cells_dict = readClonalPhylogeny(args.tree)
par2child, child2par = treeWithInternalNodes(args.tree)
print(" FILE ------->>>> ",args.tree)
print("Parent child --> ",par2child)
print("Child parent --> ",child2par)

input_df = pd.read_csv(args.genotype, sep=' ', header=None)
#print(input_df)
transposed_df = input_df.T
#print(transposed_df)
cell_genotype = transposed_df.values.tolist()[1:]
#print(cell_genotype[0])
cluster_mutation = getClonalMutation(cluster_cells_dict, cell_genotype)
updated_cluster_mutation = getNodeMutation(par2child, child2par, cluster_mutation)
mut_edge,edge_mut,parallel_edges = getEdgeMutation(par2child,updated_cluster_mutation)
sameSNVs = sameEdgeMutation(mut_edge)
#inComSNVs = parallelSNVs(edge_mut, parallel_edges)
#print(len(sameSNVs)," ",len(inComSNVs))
#print(len(set(sameSNVs.keys())-set(inComSNVs.keys())))
ancestralSNVs = getAncestralRelation(par2child, updated_cluster_mutation)
inComSNVs = getInComparableSNVs(updated_cluster_mutation, sameSNVs, ancestralSNVs)

# now think how to save these values for evaluation
save_files(set(sameSNVs.keys()),args.output+'/sameSNVs.txt')
save_files(inComSNVs,args.output+'/inComSNVs.txt')
save_files(set(ancestralSNVs.keys()),args.output+'/ancestralSNVs.txt')
