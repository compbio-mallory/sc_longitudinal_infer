import argparse
import json

#''' Read the tree file and get the SNV relations '''
#def readFile(treeFile):
#    fname = open(treeFile,'r')
#    fline = fname.readline().rstrip('\n')
#    node_mut = {}
#    while(fline != ""):
#        if "node id" in fline:
#        fline = fname.readline().rstrip('\n')

#''' save the SNV pairs to a file. '''
#def save_files(SNVset,fileName):
#    with open(fileName,'w') as f:
#        f.write(str(SNVset))

''' Find the row sum to identify the root. '''
def findRowSum(row_list):
    total = 0
    for i in row_list:
        total = total + i
    return total

''' Compare rows to connect nodes with additional 1 mutation. '''
def compareRows(B_matrix, rowIndex, colIndex, connected_nodes):
    #print(" compare rows j ",colIndex)
    nodes_to_connect = []
    for i in range(rowIndex,len(B_matrix)):
        count1s = 0 # keep a count of additional 1s. We need just one additional 1 to connect the nodes.
        lastMutIndex = 0
        for j in range(colIndex,len(B_matrix[i])):
            if B_matrix[i][j] == 1:
                count1s = count1s+1
                lastMutIndex = j
        if i not in connected_nodes and count1s == 1:
            connected_nodes.append(i)
            nodes_to_connect.append(i)

    #print(nodes_to_connect)
    return nodes_to_connect

''' Save the node mutations. Then update the dictionary to include only new mutations on the edges. '''
def getNodeMutations(B_matrix, parent_child):
    node_mut = {}
    for par, child in parent_child.items():
        #if par in node_mut:
        #    continue
        # First get the parent mutation
        par_mut_arr = B_matrix[par]
        par_mutations = []
        for j in range(len(par_mut_arr)):
            if par_mut_arr[j] == 1:
                par_mutations.append(j)

        if par not in node_mut:
            node_mut[par] = par_mutations
        # Then get the child mutations
        for c in child:
            c_mut_arr = B_matrix[c]
            child_mutations = []
            for j in range(len(c_mut_arr)):
                if c_mut_arr[j] == 1:
                    child_mutations.append(j)
            new_child_mut = set(child_mutations) - set(par_mutations)
            node_mut[c] = list(new_child_mut)

    #for par, child in parent_child.items():
    #    par_mut = node_mut[par]
    #    for c in child:
    #        c_mut = node_mut[c]
    #        new_mut = set(c_mut) - set(par_mut)
    #        node_mut[c] = new_mut
    return node_mut

''' Build the clonal tree from the perfect phylogeny matrix B. '''
def buildClonalTree(Bjson):
    #B_matrix = [[1,0,0,0,0],[1,1,0,0,0],[1,1,1,0,0],[1,1,0,1,0],[1,1,1,0,1],[1,0,1,0,0]]
    #B_matrix = [[1,0,0,0,0,0],[1,1,0,0,0,0],[1,1,1,0,0,0],[1,1,0,1,0,0],[1,0,0,0,1,0],[1,0,0,0,1,1]]
    #node_mut = {} # keep track of mutations for each nodes
    f = open(Bjson,'r')
    B_dict = json.load(f)
    B_matrix = list(B_dict.values())[0]

    parent_child = {} # Parent is the key and children are its values
    connected_nodes = [] # Use this list to exclude connected nodes. We only have 1 edge between the nodes.
    lastMutIndex = 0 # use this to update the col index or mutations
    #print(B_matrix)
    rootNode = ""

    for i in range(len(B_matrix)):
        mut_list = [] # keep track of node mutations
        row_sum = 0
        if rootNode != "":
            # First identify the root which is node 0 here
            row_sum = findRowSum(B_matrix[i])
        #print(" Connect node ",i)
        for j in range(lastMutIndex,len(B_matrix[i])):
            #print(" j index ",j)
            if row_sum == 1: # Find first node connected to the root
                rootNode = i
                child_nodes = compareRows(B_matrix, i+1, j+1, connected_nodes)
                lastMutIndex = j+1
                parent_child[i] = child_nodes
                #if B_matrix[i+1][j+1] == 1:
                #    connected_nodes.append(str(i)+'_'+str(i+1))# Connect the nodes
                break

            child_nodes = compareRows(B_matrix, i+1, j+1, connected_nodes)
            lastMutIndex = j+1
            parent_child[i] = child_nodes
            break

    node_mut = getNodeMutations(B_matrix, parent_child)
    return parent_child, node_mut

''' If there is just one branch then merge those children in one clone and update the node mutations. '''
def mergeOneBranches(parent_child, node_mut):
    merge_children = {}
    skip_nodes = []
    for par, child in parent_child.items():
        print(" Parent ",par)
        if len(child) == 1: # only check for single branches which means that there is only 1 child.
            for mp, mc in merge_children.items():
                if par in mc:
                    mc.append(child[0])
                    merge_children[mp] = mc
                    skip_nodes.append(par)
            if par not in skip_nodes:
                #print(" Parent key ",par)
                merge_children[par] = child
    print(" Merged children ",merge_children)
    print(" Skip nodes ",skip_nodes)

    updated_node_mut = {}
    updated_par_child = {}
    exclude_nodes = []
    for node, mut in node_mut.items():
        if node in merge_children:
            p_mut = node_mut[node]
            for child in merge_children[node]:
                if len(parent_child[child]) == 1:
                    exclude_nodes.append(child) # if a child has more children then don't exclude it
                p_mut.extend(node_mut[child])
                updated_node_mut[node] = list(set(p_mut))
            updated_par_child[node] = merge_children[node] 
        elif node not in exclude_nodes:
            updated_node_mut[node] = mut
        if node not in skip_nodes:
            updated_par_child[node] = parent_child[node]

    print(" Exclude nodes ",exclude_nodes)
    print(" Updated mutations ",updated_node_mut)
    print(" Updated parent child ",updated_par_child)
    
    return updated_par_child, updated_node_mut

def getNodeMutations_ForAncestral(updated_par_child, updated_node_mut):
    node_mut_forAncestral = {}
    for par, child in updated_par_child.items():
        p_mut = updated_node_mut[par]

        if par not in node_mut_forAncestral:
            node_mut_forAncestral[par] = p_mut

        for c in child:
            if c not in updated_node_mut:
                continue

            c_mut = updated_node_mut[c]
            if par in node_mut_forAncestral:
                p_mut = node_mut_forAncestral[par]

            parent_mut = list(set(p_mut))
            child_mut = list(set(c_mut))
            merged_mut = list(set(parent_mut + child_mut))
            node_mut_forAncestral[c] = merged_mut
    print(" Ancestral mutations ",node_mut_forAncestral)
    return node_mut_forAncestral

''' if pair of SNVs already checked then prevent checking. '''
def checkPairedSNVs(SNV_pairs,mut1,mut2):
    for mut in SNV_pairs:
        mut_list = mut.split('_')
        #print(mut1," ",mut2," ",mut_list)
        if str(mut1) in mut_list and str(mut2) in mut_list:
            return True

''' List SNVpairs that appear on the same branch. '''
def findSNVpairs_OnSameBranches(node_mut):
    SNV_pairs = []
    for node, mut in node_mut.items():
        for i in range(len(mut)):
            for j in range(i+1,len(mut)):
                if mut[i] == mut[j] or checkPairedSNVs(SNV_pairs, mut[i], mut[j]):
                    continue
                SNV_pairs.append(str(mut[i])+'_'+str(mut[j]))
    print("SNV pair on same branch ",set(SNV_pairs))
    return set(SNV_pairs)

''' SNVs appearing on ancestral branch. '''
def findSNVpairs_OnAncestralBranches(pID_cID, node_mut, sameSNVs):
    SNV_pairs = []
    for pid,cid in pID_cID.items():
        #if '1_-1' in pid or '1.5_0' in pid:
        #    continue
        #parent_nodeId = pid.split('_')[1]
        p_mut = node_mut[pid]
        #print(" Parent mut ",p_mut)
        for child in cid:
            if child not in node_mut:
                continue
            child_mut = node_mut[child]
            #print(" Child mut ",child_mut)
            for pmut in p_mut:
                for cmut in child_mut:
                    if pmut == cmut or checkPairedSNVs(SNV_pairs, pmut, cmut) or str(pmut)+'_'+str(cmut) in sameSNVs or str(cmut)+'_'+str(pmut) in sameSNVs:
                        continue
                    SNV_pairs.append(str(pmut)+'_'+str(cmut))
    print(" Ancestral SNVs ",set(SNV_pairs))
    return set(SNV_pairs)

''' Find incomparable SNV pairs. '''
def findSNVpairs_incomparable(node_mut, SNV_pairs_onSameBranch, SNV_pairs_OnAncestralBranch):
    SNV_pairs = []
    same_SNVpairs = list(SNV_pairs_onSameBranch)
    ancestral_SNVpairs = list(SNV_pairs_OnAncestralBranch)
    for node1, mut1 in node_mut.items():
        for node2, mut2 in node_mut.items():
            for m1 in mut1:
                for m2 in mut2:
                    if str(m1) == str(m2) or checkPairedSNVs(SNV_pairs,m1,m2) or (str(m1)+'_'+str(m2) in same_SNVpairs) or (str(m1)+'_'+str(m2) in ancestral_SNVpairs) or (str(m2)+'_'+str(m1) in same_SNVpairs) or (str(m2)+'_'+str(m1) in ancestral_SNVpairs):
                        continue
                    SNV_pairs.append(str(m1)+'_'+str(m2))
    return set(SNV_pairs)

''' save the SNV pairs to a file. '''
def save_files(SNVset,fileName):
    with open(fileName,'w') as f:
        f.write(str(SNVset))

parser = argparse.ArgumentParser()
parser.add_argument("-B", "--B",dest = "B", help="LACE's perfect phylogeny matrix.")
parser.add_argument("-output", "--output",dest = "output", help="Output path to save the files.")
args = parser.parse_args()
#readFile(args.tree)
parent_child, node_mut = buildClonalTree(args.B)
print(" Parent child ")
print(parent_child)
print(" Node mut ")
print(node_mut)
updated_par_child, updated_node_mut = mergeOneBranches(parent_child, node_mut)
node_mut_forAncestral = getNodeMutations_ForAncestral(updated_par_child, updated_node_mut)

SNVpairs_OnSameBranches = findSNVpairs_OnSameBranches(updated_node_mut)
SNVpairs_OnAncestralBranch = findSNVpairs_OnAncestralBranches(updated_par_child, node_mut_forAncestral, SNVpairs_OnSameBranches)
SNVpairs_inComp = findSNVpairs_incomparable(updated_node_mut, SNVpairs_OnSameBranches, SNVpairs_OnAncestralBranch)
# Save the SNVpairs files.
save_files(SNVpairs_OnSameBranches,args.output+'/sameSNVs.txt')
save_files(SNVpairs_OnAncestralBranch,args.output+'/ancestralSNVs.txt')
save_files(SNVpairs_inComp,args.output+'/inComSNVs.txt')
