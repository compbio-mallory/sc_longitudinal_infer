import argparse

''' Get parent child nodes here. '''
def readTree(tree):
    tfile = open(tree,'r')
    tline = tfile.readline().rstrip('\n') # Skipping the 1st line
    tline = tfile.readline().rstrip('\n') # Skipping the 2nd line
    tline = tfile.readline().rstrip('\n')
    par_child = {}
    child_list = []
    while(tline != ""):
        if "node [color=lightgrey" in tline:
            break
        par_child_arr = tline.split(' -> ')
        par = str(int(par_child_arr[0]) - 1)
        child = str(int(par_child_arr[1].rstrip(';')) - 1)
        child_list.append(child)
        if par in par_child:
            t_child = par_child[par]
            t_child.append(child)
            par_child[par] = t_child
        else:
            t_child = []
            t_child.append(child)
            par_child[par] = t_child
        tline = tfile.readline().rstrip('\n')
    print(par_child)
    print(" Children list ",child_list)
    return par_child, child_list

''' Arrange the nodes on the tree. '''
def arrangeTree(par_child, child_list):
    clonal_tree = {} # Have parent clones as key and child clones as values
    clone_mut = {} # Have clones as key and mutations as values
    merge_subclones = {} # Record subclones that need to be merged
    par_list = list(par_child.keys())

    clone_number = 0
    merge_children = []
    for par, child in par_child.items():
        #mut_list = []
        clone_list = []
        par_clone_number = 'C'+par

        if par not in child_list: # Get the root node
            print(" Root is ",par)
            clone_list.append(par_clone_number)
            clonal_tree['Root'] = clone_list
            clone_mut['Root'] = [] # Root doesn't have a mutations so this is empty
            #mut_list.extend(child)
            clone_mut[par_clone_number] = [par]
            #continue

        if len(child) == 1:
            #print(" Merge children ",merge_children)
            temp_list = []
            temp_list.append(par+'_'+child[0])
            merge_children.append(temp_list)

        for c in child:
            if c not in par_list: # Merge the multiple subclones here to a single parent
                if par in merge_subclones:
                    t_child = merge_subclones[par]
                    t_child.append(c)
                    merge_subclones[par] = t_child
                else:
                    t_child = []
                    t_child.append(c)
                    merge_subclones[par] = t_child
            else:
                child_clone_number = 'C'+c
                #mut_list.append(c)
                if par_clone_number in clonal_tree:
                    t_clone = clonal_tree[par_clone_number]
                    t_clone.append(child_clone_number)
                    clonal_tree[par_clone_number] = t_clone
                    clone_mut[child_clone_number] = [c] # Each child exists as an individual clone if not merged
                else:
                    t_clone = []
                    t_clone.append(child_clone_number)
                    clonal_tree[par_clone_number] = t_clone
                    clone_mut[child_clone_number] = [c] # Each child exists as a individual clone if not merged

        if par in merge_subclones:
            #mut_list.extend(merge_subclones[par])
            clone_mut[par_clone_number] = [par] # First add the parent mutations 
            child_mut_list = merge_subclones[par] # Then add the child mutations
            clone_mut['MC'+str(clone_number)] = child_mut_list
            if par_clone_number in clonal_tree:
                temp = clonal_tree[par_clone_number]
                temp.append('MC'+str(clone_number))
                clonal_tree[par_clone_number] = temp
            else:
                clonal_tree[par_clone_number] = ['MC'+str(clone_number)] # Then connect parent clone to merged child clone

        clone_number = clone_number + 1
            
    print(" Merged subclones ",merge_subclones)
    print(" Clone mut ",clone_mut)
    print(" Clone Tree ",clonal_tree)
    print(" Merged children ",merge_children)
    return clone_mut, clonal_tree

''' Since the tree is not ordered we just run the ancestral mutation loop second time to get the appropriate SNVs. '''
def getAncestralSNVs(clonal_tree, node_mut_forAncestral):
    for par, child in clonal_tree.items():
        if par in node_mut_forAncestral:
            p_mut = node_mut_forAncestral[par]

        parent_mut = list(set(p_mut))
        for c in child:
            c_mut = node_mut_forAncestral[c]
            child_mut = list(set(c_mut))
            merged_mut = list(set(parent_mut + child_mut))
            print(c," Merged mut ",merged_mut)
            node_mut_forAncestral[c] = merged_mut

    return node_mut_forAncestral

''' Get the edges and edge mutation from the tree. '''
def getEdges(clone_mut, clonal_tree):
    edge_mut = {}
    parallel_edges = []
    node_mut_forAncestral = {} # contains node mutations as it appears

    for par, child in clonal_tree.items(): 
        par_mut = clone_mut[par]
        #print(" Parent mut ",par_mut)
        if par not in node_mut_forAncestral:
            node_mut_forAncestral[par] = par_mut
        
        #print(" Parent ",par," Child ",child)
        node_parallel_edges = []
        for c in child:
            c_mut = clone_mut[c]
            #print(" Child mut ",c_mut)
            edge = par+'_'+c
            edge_mut[edge] = list(set(c_mut) - set(par_mut))
            node_parallel_edges.append(edge)

            if par in node_mut_forAncestral: # get updated parent mutations
                p_mut = node_mut_forAncestral[par]
                
            parent_mut = list(set(p_mut))
            child_mut = list(set(c_mut))
            merged_mut = list(set(parent_mut + child_mut))
            #print(c," Merged mut ",merged_mut)
            node_mut_forAncestral[c] = merged_mut

        if len(node_parallel_edges) > 1:
            parallel_edges.append(node_parallel_edges)
    print(" Edge mutations ")
    print(edge_mut)
    print(" Parallel mutations ",parallel_edges)
    final_node_mut_forAncestral = getAncestralSNVs(clonal_tree, node_mut_forAncestral)
    #print(" Ancestral Node mutations ",node_mut_forAncestral)

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
    return mut_edge, edge_mut, parallel_edges, final_node_mut_forAncestral
        
''' If there is just one branch then merge those children in one clone and update the node mutations. '''
def mergeOneBranches(parent_child, node_mut):
    rootNode = parent_child['Root']
    rootNode_children = parent_child[rootNode[0]]
    
    updated_clonal_tree = {}
    merge_children = []
    parent_child_keys = list(parent_child.keys())
    #not_checked_nodes = [] # use this list to recheck some of the not added nodes

    for rc in rootNode_children:
        nodes_to_merge = []
        not_checked_nodes = [] # use this list to recheck some of the not added nodes
        nodes_to_merge.append(rc)

        should_restart = True
        while should_restart:
            should_restart = False
            for i in range(len(parent_child_keys)):
                #print(" Value of i ",i)
                #print(" Value of rc ",rc)
                if rc == parent_child_keys[i]:
                    child = parent_child[rc]
                    if len(child) == 1:
                        if child[0] in parent_child:
                            child_child = parent_child[child[0]] # check if a node has more children or merged subclone. Then skip them from adding.
                            if 'MC' in child_child[0] or len(child_child) > 1:
                                print(" Not added node ",child)
                                #not_checked_nodes.append(child)
                                if 'MC' not in child_child[0]:
                                    rootNode_children.append(child_child[0]) # Since this branches out make it a root
                                break
                            #if len(child_child) > 1:
                            #    for cc in child_child: # There might be inner nodes
                            #        inner_nodes_to_merge = []
                            #        if cc == 'MC':
                            #            print(" Not added MC node ",cc)
                            #        else:
                            #            inner_nodes_to_merge.append(cc)
                            #            rc = cc
                            #            should_restart = True # start the loop from the beginning
                            #            break
                        nodes_to_merge.append(child[0])
                        rc = child[0]
                        should_restart = True # start the loop from the beginning
                        break
        merge_children.append(nodes_to_merge)
    print(" Updated merge nodes ",merge_children)
    #print(" Inner nodes to merge ",inner_nodes_to_merge)

    for i in range(len(merge_children)): 
        firstNode = merge_children[i][0] # to serve as the parent
        lastNode = merge_children[i][len(merge_children[i])-1] # to connect the nodes
        name_mergedNodes = 'M'+lastNode

        parent_child[firstNode] = [name_mergedNodes]
        if lastNode in parent_child:
            parent_child[name_mergedNodes] = parent_child[lastNode] # Connect the nodes 

        node_mut[firstNode] = [k.split('C')[1] for k in merge_children[i][:]]
        node_mut[name_mergedNodes] = node_mut[lastNode]

        for j in range(1,len(merge_children[i])):
            if merge_children[i][j] in parent_child:
                parent_child.pop(merge_children[i][j]) # Remove merged children from parent_child
            node_mut.pop(merge_children[i][j]) # Remove their mutations as well

    print(" Updated clonal tree ",parent_child) 
    print(" Updated node mutation ",node_mut)

    return parent_child, node_mut

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
def getAncestralRelation(par2child, updated_cluster_mutation, sameSNVs):
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
                    if checkPairedSNVs(ancestralSNVs,m1,m2) or m1 == m2 or str(m1)+'_'+str(m2) in sameSNVs or str(m2)+'_'+str(m1) in sameSNVs:
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
    print(" Ancestral SNVs ",ancestralSNVs.keys())
    return ancestralSNVs

''' Get the SNV pairs that are incomparable. '''
def getInComparableSNVs(clone_mut, sameSNVs, ancestralSNVs):
    inCompSNVs = []
    sameSNVs_list = set(sameSNVs.keys())
    ancestralSNVs_list = set(ancestralSNVs.keys())
    print(" sameSNVs no. ",len(sameSNVs_list)," ancestralSNVs no. ",len(ancestralSNVs_list))
    for node1, mut1 in clone_mut.items():
        for node2, mut2 in clone_mut.items():
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
        
parser = argparse.ArgumentParser()
parser.add_argument("-tree", "--tree",dest ="tree", help="SCITE's inferred tree GV.")
parser.add_argument("-output", "--output",dest ="output", help="Output dir to save the files.")
args = parser.parse_args()
par_child, child_list = readTree(args.tree)
clone_mut, clonal_tree = arrangeTree(par_child, child_list)
updated_clonal_tree, updated_clone_mut = mergeOneBranches(clonal_tree, clone_mut)
mut_edge, edge_mut, parallel_edges, node_mut_forAncestral = getEdges(updated_clone_mut, updated_clonal_tree)

#print(" Merge branches code ")
#mergeOneBranches(clonal_tree, clone_mut)
print("Mutation on edges ",edge_mut)
sameSNVs = sameEdgeMutation(mut_edge)
#inComSNVs = parallelSNVs(edge_mut, parallel_edges)
ancestralSNVs = getAncestralRelation(updated_clonal_tree, node_mut_forAncestral, set(sameSNVs.keys()))
inComSNVs = getInComparableSNVs(clone_mut, sameSNVs, ancestralSNVs)

save_files(set(sameSNVs.keys()),args.output+'/sameSNVs.txt')
#save_files(set(inComSNVs.keys()),args.output+'/inComSNVs.txt')
save_files(inComSNVs,args.output+'/inComSNVs.txt')
save_files(set(ancestralSNVs.keys()),args.output+'/ancestralSNVs.txt')
