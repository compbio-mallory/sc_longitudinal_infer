import argparse
import json
import copy

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
            if int(B_matrix[i][j]) == 1:
                count1s = count1s+1
                lastMutIndex = j
        if i not in connected_nodes and count1s == 1:
            connected_nodes.append(i)
            nodes_to_connect.append(i)

    #print(nodes_to_connect)
    return nodes_to_connect

''' Build the clonal tree from the perfect phylogeny matrix B. '''
def buildClonalTree(B_matrix, snv_array):
    #B_matrix = [[1,0,0,0,0],[1,1,0,0,0],[1,1,1,0,0],[1,1,0,1,0],[1,1,1,0,1],[1,0,1,0,0]]
    #B_matrix = [[1,0,0,0,0,0],[1,1,0,0,0,0],[1,1,1,0,0,0],[1,1,0,1,0,0],[1,0,0,0,1,0],[1,0,0,0,1,1]]
    #node_mut = {} # keep track of mutations for each nodes
    #f = open(Bjson,'r')
    #B_dict = json.load(f)
    #B_matrix = list(B_dict.values())[0]

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
            #if i >= len(B_matrix[i])-1 or j >= len(B_matrix[i])-1:
            #    continue
            #print(i," j index ",j," i+1 ",i+1," j+1 ",j+1)
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

     # To avoid confusion also update the parent_child to have values from the snv array
    new_parent_child = {}
    for p, child in parent_child.items():
        parent_from_snvs = snv_array[p]
        c_list = []
        for c in child:
            c_list.append(snv_array[c])
        new_parent_child[parent_from_snvs] = c_list
    print("New parent child list ",new_parent_child)

    #node_mut = getNodeMutations(B_matrix, parent_child)
    #new_node_mut = {}
    ## Update the node mutation array to have original mutations
    #for n, mut in node_mut.items():
    #    mut_list = []
    #    for m in mut:
    #        mut_list.append(snv_array[m])
    #    new_node_mut[n] = mut_list

    #print("New node mutation ",new_node_mut)
    return new_parent_child

''' Build the clonal tree from the perfect phylogeny matrix B. '''
#def buildClonalTree(Bmatrix, snv_array):
#    parent_child = {} # Parent is the key and children are its values
#    connected_nodes = []
    # First row is Root so we ignore it
#    for i in range(1,len(Bmatrix)-1):
        #child_nodes = compareRows(Bmatrix, i+1, i+1, [])
#        child_nodes = []
#        for j in range(len(Bmatrix[i])):
#            nodes = compareRows(Bmatrix, i+1, j+1, [])
#            child_nodes.extend(nodes)
            #print(child_nodes)
            #if int(Bmatrix[i+1][j+1]) == 1:
            #    child_nodes.append(snv_array[i+1])
#        print(child_nodes)

#        parent_child[snv_array[i]] = child_nodes
#    print(parent_child)

''' Read the csv file for perfect phylogenetic matrix B and return the B_matrix where first row represents the Root and SNVs belonging to each array index of B.. '''
def readBfile(Bcsv):
    f = open(Bcsv,'r')
    fline = f.readline().rstrip('\n') 
    snv_array = {} # Create a dict with array indices -> SNVs
    B_matrix = []
    header = True
    while(fline != ""):
        if header:
            larr = fline.split(",")[2:]
            print(larr)
            for i in range(len(larr)):
                tempArr = larr[i].split('V')
                if len(tempArr) > 2:
                    snv_array[i] = int(tempArr[len(tempArr)-1].rstrip('"'))
                else:
                    snv_array[i] = int((larr[i].split('V')[1]).rstrip('"'))
            header = False
            fline = f.readline().rstrip('\n')
            continue
        
        if "Root" in fline.split(",")[0]: # No need to add the Root for evaluation
            print("Encountered Root")
            fline = f.readline().rstrip('\n')
            continue

        larr = fline.split(",")[2:]
        larr = [int(i) for i in larr]
        B_matrix.append(larr)
        fline = f.readline().rstrip('\n')
    print("Array indices to SNVs ",snv_array)
    print("B_matrix ",B_matrix)
    print(len(B_matrix)," ",len(B_matrix[0]))
    total_mutations = len(B_matrix)
    return B_matrix, snv_array, total_mutations

''' If there is just one branch then merge those children in one clone and update the node mutations. '''
def mergeOneBranches(parent_child):
    parent_child_dict = copy.deepcopy(parent_child)
    merge_children = {}
    skip_nodes = []
    for par, child in parent_child_dict.items():
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
    
    return merge_children

''' if pair of SNVs already checked then prevent checking. '''
def checkPairedSNVs(SNV_pairs,mut1,mut2):
    for mut in SNV_pairs:
        mut_list = mut.split('_')
        #print(mut1," ",mut2," ",mut_list)
        if str(mut1) in mut_list and str(mut2) in mut_list:
            return True

''' List SNVpairs that appear on the same branch. '''
# Input node_mut = merged mutations based on parent child relationship
def findSNVpairs_OnSameBranches(node_mut):
    SNV_pairs = []
    for node, mut in node_mut.items():
        for i in range(len(mut)):
            for j in range(i+1,len(mut)):
                if mut[i] == mut[j] or checkPairedSNVs(SNV_pairs, mut[i], mut[j]):
                    continue
                SNV_pairs.append(str(mut[i])+'_'+str(mut[j]))
            SNV_pairs.append(str(node)+'_'+str(mut[i]))
    print("SNV pair on same branch ",SNV_pairs)
    return set(SNV_pairs)

''' SNVs appearing on ancestral branch. '''
#def findSNVpairs_OnAncestralBranches(pID_cID, node_mut, sameSNVs):
#    SNV_pairs = []
#    for pid,cid in pID_cID.items():
        #if '1_-1' in pid or '1.5_0' in pid:
        #    continue
        #parent_nodeId = pid.split('_')[1]
#        p_mut = node_mut[pid]
        #print(" Parent mut ",p_mut)
#        for child in cid:
#            if child not in node_mut:
#                continue
#            child_mut = node_mut[child]
            #print(" Child mut ",child_mut)
#            for pmut in p_mut:
#                for cmut in child_mut:
#                    if pmut == cmut or checkPairedSNVs(SNV_pairs, pmut, cmut) or str(pmut)+'_'+str(cmut) in sameSNVs or str(cmut)+'_'+str(pmut) in sameSNVs:
#                        continue
#                    SNV_pairs.append(str(pmut)+'_'+str(cmut))
#    print(" Ancestral SNVs ",set(SNV_pairs))
#    return set(SNV_pairs)

# Input is to check if current parent is a child in merged_mutations to get all the mutations for ancestral relationship
# Output is return the parent and child mutations to build all the SNV pairs
def getAllParentMutations(p, parent_child, merged_mutations):
    all_mut = []
    for parent, child in parent_child.items():
        if p in child:
            # Then look for the mutations for this parent in merged_mutations to get all the ancestral mutations
            for pm, cm in merged_mutations.items():
                if parent in cm:
                    all_mut.append(pm)
                    all_mut.extend(cm)
    return all_mut

def findSNVpairs_OnAncestralBranches(parent_child, merged_mutations):
    # {5: [2], 2: [6], 6: [0], 0: [1], 1: [3], 3: [4], 4: [12, 9], 12: [11], 11: [13], 13: [16, 7], 16: [14], 14: [17], 17: [15], 15: [18], 18: [19], 19: [], 9: [8], 8: [10], 10: [], 7: []}
    # {5: [2, 6, 0, 1, 3, 4], 12: [11, 13], 16: [14, 17, 15, 18, 19], 9: [8, 10]}
    # Is 5 from above child of anyone's in parent_child? No then great do nothing
    # If 12 is child of anyone from the parent_child then note its parent then find all mutations associated with it. (Already implemnted this)
    # Update a dictionary with 12, 11, 13 as keys and values as the mutations to be considered for ancestral relationships. 
    # Now for node 16 first check if 13 exists in the above dictionary. If yes then simply add those mutations for ancestral relationships and update above dictionary with '16, 14, 15, ..' as keys.
    # Ex of the dictionary: {'12,11,13': '0,1,2,3,4,5,6', '16,14,17,15,18,19': '11,12,13,0,1,2,3,4,5,6', '9,8,10': '0,1,2,3,4,5,6' 
    # Now compare the keys with the values to get the ancestral SNVs
    SNV_pairs = []
    ancestral_mut = {}
    for parent, child in merged_mutations.items():
        #print("Parent ",parent," child ",child)
        allMut = getAllParentMutations(parent, parent_child, merged_mutations)
        if allMut != []:
            #print("More mutations ",allMut)
            # Check if there are existing parent in the ancestral_mut then add their parents too
            if ','.join(map(str,allMut)) in ancestral_mut:
                allMut.extend(ancestral_mut[','.join(map(str,allMut))])

            key = str(parent)+','+','.join(map(str,child))
            ancestral_mut[key] = allMut

    # Some mutations like 7 (node with single mutation and no children) might not be in merged mutations so we need to make sure it is also added in the ancestral_mut.
    # We check the original parent_child to find those mutations first
    #print("Parent child ",parent_child)
    for p, clist in parent_child.items():
        if len(clist) > 1: # single all single branches are merged we don't need to check any further
            for c in clist:
                if c not in merged_mutations:
                    #print("Missed mutation ",c)
                    allMut = getAllParentMutations(c, parent_child, merged_mutations)
                    if allMut != []:
                        # Check if there are existing parent in the ancestral_mut then add their parents too
                        if ','.join(map(str,allMut)) in ancestral_mut:
                            allMut.extend(ancestral_mut[','.join(map(str,allMut))])
                        ancestral_mut[str(c)] = allMut

    print("Ancestral mut dict ",ancestral_mut)
    for k, v in ancestral_mut.items():
        k_list = k.split(',')
        for mut1 in k_list:
            for mut2 in v:
                SNV_pairs.append(str(mut1)+'_'+str(mut2))
    print("Ancestral SNV pairs ",SNV_pairs)
    return set(SNV_pairs)

''' Get all possible combinations of SNVs. '''
def allSNVpairs(totalMut):
    SNV_pairs = []
    for i in range(totalMut):
        for j in range(1, totalMut):
            if i == j or (str(i)+'_'+str(j) in SNV_pairs) or (str(j)+'_'+str(i) in SNV_pairs):
                continue
            SNV_pairs.append(str(i)+'_'+str(j))

    print("All possible SNVs ",SNV_pairs)
    return SNV_pairs

''' Find incomparable SNV pairs. '''
def findSNVpairs_incomparable(allSNVs, SNV_pairs_onSameBranch, SNV_pairs_OnAncestralBranch):
    SNV_pairs = copy.deepcopy(allSNVs)
    same_SNVpairs = list(SNV_pairs_onSameBranch)
    ancestral_SNVpairs = list(SNV_pairs_OnAncestralBranch)
    for snvs in allSNVs:
        m1 = snvs.split('_')[0]
        m2 = snvs.split('_')[1]
        if (str(m1)+'_'+str(m2) in same_SNVpairs) or (str(m1)+'_'+str(m2) in ancestral_SNVpairs) or (str(m2)+'_'+str(m1) in same_SNVpairs) or (str(m2)+'_'+str(m1) in ancestral_SNVpairs):
            SNV_pairs.remove(snvs)

    print("Incomparable SNVs ",SNV_pairs)
    return set(SNV_pairs)

''' save the SNV pairs to a file. '''
def save_files(SNVset,fileName):
    with open(fileName,'w') as f:
        f.write(str(SNVset))

parser = argparse.ArgumentParser()
parser.add_argument("-B", "--B",dest = "B", help="LACE's perfect phylogeny matrix.")
parser.add_argument("-output", "--output",dest = "output", help="Output path to save the files.")
args = parser.parse_args()
Bmatrix, snv_array, totalMut = readBfile(args.B)
#buildClonalTree(Bmatrix, snv_array)
#readFile(args.tree)
parent_child = buildClonalTree(Bmatrix, snv_array)
print(" After buildClonalTree Parent child ")
print(parent_child)
#print(" Node mut ")
#print(node_mut)
merged_mutations = mergeOneBranches(parent_child)
print("After merged mutations ",parent_child)
allSNVs = allSNVpairs(totalMut)
#node_mut_forAncestral = getNodeMutations_ForAncestral(updated_par_child, updated_node_mut)

SNVpairs_OnSameBranches = findSNVpairs_OnSameBranches(merged_mutations)
#SNVpairs_OnSameBranches = []
SNVpairs_OnAncestralBranch = findSNVpairs_OnAncestralBranches(parent_child, merged_mutations)
SNVpairs_inComp = findSNVpairs_incomparable(allSNVs, SNVpairs_OnSameBranches, SNVpairs_OnAncestralBranch)
# Save the SNVpairs files.
save_files(SNVpairs_OnSameBranches,args.output+'/sameSNVs.txt')
save_files(SNVpairs_OnAncestralBranch,args.output+'/ancestralSNVs.txt')
save_files(SNVpairs_inComp,args.output+'/inComSNVs.txt')
