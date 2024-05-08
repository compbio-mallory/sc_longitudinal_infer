import argparse
import copy

''' Define the tree structure. '''
class Node:
    def __init__(self, id_):
        self.id = id_  # node or child id
        self.pID = -1 # parent node
        self.mutations = [] # mutations incoming on this node
        self.edges = "" # incoming edges
        self.children = [] # child nodes
        self.timepoint = -1
        self.cells = [] # save the cells belonging to this node
        self.type = "root" # We start with the root node

''' Print the contents of the Tree after updating. '''
def printTree(Tree):
    print("=============== TREE =====================")
    for j in range(len(Tree)):
        print(" ID ",Tree[j].id," PID ",Tree[j].pID," Timepoint ",Tree[j].timepoint," Mutations ",Tree[j].mutations," Type ",Tree[j].type," Edges ",Tree[j].edges," Children ",Tree[j].children," Cells ",len(Tree[j].cells))
    print(" ====================================== ")

''' Given a genotype returns an array with indices of 1s. '''
def getMutations(genotype):
    mut = []
    for i in range(len(genotype)):
        if genotype[i] == 1:
            mut.append(i)
    return mut

''' Build the graph after correcting the clusters based on likelihood values. '''
# Input is clusters and its assigned cells along with its genotype at each timepoint
# Updates the initial Tree from the clustering results
def buildTree(tp_cluster_cells, tp_cluster_genotype):
    clone_node = {} # Map the clones to the nodes in Graph for downstream analysis.
    Tree = []
    Tree.append(Node(0)) # Append the root node.
    Tree[0].id = 0
    node = 1
    tp_nodes = {} # Can be used later to check consecutive nodes at each timepoint
    tp_nodes[-1] = [0] # Include root node
    print(" Inside buildTree ============================== ")
    #print("tp_cluster_cells ",tp_cluster_cells)
    for tp, clusters in tp_cluster_cells.items():
        for cluster, cell in clusters.items():
            Tree.append(Node(node))
            Tree[node].id = node
            Tree[node].timepoint = tp
            Tree[node].cells = cell
            Tree[node].type = "subclone"
            Tree[node].mutations = getMutations(tp_cluster_genotype[tp][cluster])
            #clone_node[str(cluster)+"_"+str(tp)] = node
            clone_node[str(tp)+"_"+str(cluster)] = node
            if tp in tp_nodes:
                tp_nodes[tp].append(node)
            else:
                tp_nodes[tp] = [node]
            node = node+1
    print(" TP nodes ",tp_nodes)
    return Tree, clone_node, tp_nodes

''' Connect the nodes across timepoints having the same genotype. '''
# Input is Tree
def connectSameNodes(Tree):
    connectedNodes = []
    for i in range(len(Tree)):
        for j in range(len(Tree)):
            if i == j:
                continue
            if Tree[i].timepoint + 1 == Tree[j].timepoint:
                if Tree[i].mutations == Tree[j].mutations:
                    Tree[i].children.append(j)
                    Tree[j].pID = i
                    connectedNodes.append(j)
    return Tree, connectedNodes

''' Build the hash table with mutations as key and nodes as values to reduce complexity while listing the subsets. 
While building the hash table we don't include mutations having one nodes as they had to be eliminated later anyway. '''
# Input is all mutations of nodes in next timepoint
def build_hash_table(node_mut):
    #mut_dict = {} # Hash table 1 where mutation is the key and nodes are values
    hash_table = {} # Hash table with merged mutations as keys and nodes as values
    #added_nodes = []
    for n1, m1 in node_mut.items():
        #if n1 in added_nodes:
        #    continue
        for n2, m2 in node_mut.items():
            if n1 == n2:
                continue
            if m1 == m2:
                mut_key = ",".join(map(str,m1))
                if mut_key in hash_table:
                    hash_table[mut_key].append(n2)
                else:
                    hash_table[mut_key] = [n1,n2]
                    #added_nodes.extend([n1,n2])
            elif set(m1) & set(m2):
                subset_mut = set(m1) & set(m2)
                #print(subset_mut)
                mut_key = ",".join(map(str,subset_mut))
                #added_mut.extend(list(subset_mut))
                if mut_key in hash_table:
                    hash_table[mut_key].append(n2)
                else:
                    hash_table[mut_key] = [n1,n2]
                    #added_nodes.extend([n1,n2])

    print(" HASH TABLE ============================")
    for m, n in hash_table.items():
        hash_table[m] = list(set(n))
        print(m," --> ",set(n))

    return hash_table

# Check for unobserved subclones b/w two time points.
# Pseudocode to check for unobserved subclones.
#for t1, nodes in tp_nodes.items():
#  nodes_list = nodes
#  while(usc_node_list != []):
#    t2_nodes = tp_nodes[t1+1]
#    build_hashTable(t2_nodes_mut)
    # usc_node_list = Eliminate + select unobserved subclones(requires the nodes_list) b/w t1 and t2. Save them (use the one I thought)
    #and connect them or recheck before connecting.
    # nodes_list.append(uscs)

''' Check if any node in time t has same set of mutations as any of the unobserved subclones
and eliminate those unobserved subclones. Also, eliminate unobserved subclones having 1 child node. '''
# Input is Tree, nodes in previous timepoint to check mutations, hash table with mutation as key and nodes as values.
def eliminate_usc_sets(Tree, nodes, hash_table):
    print("====== Inside eliminating unobserved subclones ==========")
    hash_table_copy = copy.deepcopy(hash_table)
    for mut, nl in hash_table.items():
        mut_list = mut.split(',')
        mut_list = [int(i) for i in mut_list] 
        present = False
        # Check if any nodes already added as unobserved subclone between consecutive timepoints.
        # If yes then don't enroll it any further.
        #if nodes_added != []:
        #    updatedNodes = filterUnobservedSubclones(nodes_added, nl)
        #    if updatedNodes == []:
        #        present = True

        if len(nl) == 1 or nl == []: # check if there is only 1 child or no children
            present = True
        else:
            for n in nodes:
                if mut_list == Tree[n].mutations:
                    present = True

        if present:
            del hash_table_copy[mut]
        else:
            continue
    return hash_table_copy

''' After enrolling one unobserved subclone remove their children from other potential unobserved subclones. '''
# Input n1 = already added nodes, n2 = nodes present in the hash table
def filterUnobservedSubclones(n1, n2):
    for n in n1:
        if n in n2:
            n2.remove(n)
    return n2

''' Recheck the parent node. If the selected parent is already connected to any of the unobserved nodes in the list of nodes
then recheck its subset of mutations and connect the two unobserved nodes. '''
# Input is Tree, parent node, list of nodes in previous timepoint and mutations of unobserved subclone to be enrolled
def recheckParentNode(Tree, node, t1_nodes, mut):
    for c in Tree[node].children:
        if Tree[c].type == "unobserved" and c in t1_nodes:
            #print(" Unobserved subclone ",c)
            #print(Tree[c].mutations," ",mut)
            if len(Tree[c].mutations) != len(mut) and set(Tree[c].mutations).issubset(set(mut)):
                #print(" Unobserved subclone ",c)
                return c
    return ""

''' We first decreasingly sort the remaining merged keys by the number of mutations they overlap with the consensus genotype of any of the nodes in time k, 
and designated the unobserved node according to this order. The designated unobserved node is linked to the node in time k whose consensus genotype has the 
maximum overlap with the unobserved node.'''
''' From the hash table select the unobserved subclone which has the largest number of child nodes.
To connect with the parent we choose the node in time t that has largest intersection of mutations. '''
# Input is the Tree, nodes in previous timepoint, hash_table, timepoint for unobserved subclone, iter No., children if list not from hash table
def select_unobservedSubclones(Tree, nodes, sorted_hash_table, usc_timepoint, iterNo, children_list):
    noOfIntersections = 0
    parentNode = {} # used to save the parentNode. There can be only one parent.
    nodes_added = []
    unobserved_subclones = []

    for mut, nodeList in sorted_hash_table.items(): # sorted hash table is a tuple here
        # re-check if this unobserved subclone should be eliminated before proceeding.
        print("Nodes added ",nodes_added," nodes list ",nodeList)
        if nodes_added != []:
            updatedNodes = filterUnobservedSubclones(nodes_added, nodeList)
            print(" Updated nodes ",updatedNodes)
            # If there is no more children left to add then don't enroll this unobserved subclone
            if updatedNodes == []:
                continue

        mut_list = mut.split(',')
        mut_list = [int(i) for i in mut_list] # Convert the mutations to integers
        print("For mutations ",mut_list)
        for n in nodes:
            #print(" Node in prev timepoint ",n)
            if len(nodes) == 1 and n == 0: # Root node doesn't have any mutations so no need to check intersections
                parentNode[0] = n
                continue
            if set(Tree[n].mutations) == set(mut_list): # Don't enroll unobserved subclones that are same as previous node
                continue
            mutSubset = set(Tree[n].mutations) & set(mut_list)
            #print(" Mutation subset ",mutSubset," ",Tree[n].mutations," ",mut_list)
            if len(mutSubset) > noOfIntersections:
                noOfIntersections = len(mutSubset)
                print("Could have selected parent ",n," overlapped mutations ",noOfIntersections)
                if iterNo > 1: # If unobserved subclone from other iterations have less mutations than existing one then don't enroll them
                    if len(mut_list) != len(Tree[n].mutations) and set(mut_list).issuperset(set(Tree[n].mutations)): # new unobserved subclone should also be a proper superset of its parent
                        print("Selected parent ",n)
                        parentNode[0] = n # Update the node with maximum no. of intersections
                else:
                    parentNode[0] = n # Update the node with maximum no. of intersections
        
        if parentNode == {}:
            continue

        newParent = recheckParentNode(Tree, parentNode[0], nodes, mut_list) # Recheck parent for multiple iterations of unobserved subclones
        if newParent != "":
            print("New parent ",newParent)
            parentNode[0] = newParent

        nodes_added.extend(nodeList)
        print("Nodes added ",nodes_added," Node list ",nodeList)
        print("unobserved_subclones ",unobserved_subclones)
        if len(nodeList) == 1: # prevent enrolling unobserved subclones having 1 child
            continue
        
        # Check before enrolling the subclones
        #checkEnroll = checkToEnroll(Tree, mut_list, parentNode[0], nodeList)
        #if checkEnroll == False:
        #    continue

        # Update the tree with parent ID and unobserved subclone and get the enrolled nodeId
        newNodeId = enrollUnobservedSubclones(Tree, mut_list, parentNode[0], nodeList, usc_timepoint, children_list)
        unobserved_subclones.append(newNodeId)
        print(" Unobserved subclone ",mut," children ",nodeList," parents ",parentNode[0]," ",Tree[parentNode[0]].children)
        print(" New node Id ",newNodeId)
    return unobserved_subclones

''' While enrolling subclone check if its children previously added to any unobserved subclones. 
If yes then check if that was its parent too. If not then don’t make this unobserved subclones its parent.
Then check if this unobserved subclone will have only one child or not child and hence shouldn’t be added.'''
#def checkToEnroll(Tree, mut, parentNode, cIDs):
#    enroll = True
#    removeChild = []
#    for c in cIDs:
#        pID = Tree[c].pID
#        if parentNode != pID and Tree[pID].type == "unobserved":
            # This child is already added so no need to add it again
#            removeChild.append(c)

#    for rc in removeChild:
#        cIDs.remove(rc)

#    if cIDs == [] or len(cIDs) == 1:
#        enroll = False
#    return enroll

''' Update the Tree to include the unobserved subclones. '''
# Input is Tree, mutations of the unobserved subclone, parent IDs and children
def enrollUnobservedSubclones(Tree, mut, pID, cIDs, usc_timepoint, children_list):
    if children_list != []:
        print("Inside Enroll unobserved subclones ")
        cIDs = children_list
    print(cIDs)
    node = len(Tree)
    Tree.append(Node(node))
    Tree[node].id = node
    Tree[node].timepoint = usc_timepoint
    Tree[node].mutations = mut
    Tree[node].pID = pID
    Tree[node].children = cIDs
    Tree[node].type = "unobserved"
    Tree[pID].children.append(node)
    print("Unobserved subclone ",node," children ",cIDs)
    # If parent has the same set of children then remove them
    for c in cIDs:
        if Tree[c].timepoint != usc_timepoint:
            print("Timepoints different !!! ")
            Tree[node].children.remove(c)
            continue
        if c in Tree[pID].children:
            #print(" pID ",pID," c ",c)
            Tree[pID].children.remove(c)
    #if pID == 0: # For root node
    #    Tree[pID].children.append(node)
    # Connect the children to its parent
    for c in cIDs:
        # Before connecting children to its parents check if there are existing parent
        # Remove the children from those previous pIDs before enrolling new children
        if Tree[c].pID != -1 and c in Tree[Tree[c].pID].children:
            Tree[Tree[c].pID].children.remove(c)
        if Tree[Tree[c].pID].children == []:
            Tree[Tree[c].pID].pID = -2
            Tree[Tree[c].pID].edges = ""
        Tree[c].pID = node
    return node

''' Gets a dictionary with node as key and mutations as values for a timepoint. '''
# Input is Tree, nodes in previous timepoints, already connected nodes with same set of mutations
def nodeMutations(Tree, nodes, connectedNodes):
    node_mut = {}
    for n in nodes:
        if n in connectedNodes: # Skip adding nodes that are already connected with same set of mutations in prev timepoint
            continue
        node_mut[n] = Tree[n].mutations
    return node_mut

''' Check if the updated hash table has tie in terms of children. 
If yes then re-sort the hash table according to no. of mutations. '''
def checkHashTable(hash_table):
    vlen = set()
    for k, v in hash_table.items():
        vlen.add(len(v))
    if len(vlen) == 1: # this would mean all mutations have same no. of children
        mut_list = list(hash_table.keys())
        mut_dict = {}
        for i in range(len(mut_list)): # convert the str to list elements to check len
            mut_len = mut_list[i].split(',')
            mut_dict[mut_list[i]] = mut_len
        sorted_mut_dict = dict(sorted(mut_dict.items(), key=lambda items:len(items[1]), reverse=True))
        new_hash_table = {}
        for m in sorted_mut_dict:
            new_hash_table[m] = hash_table[m]
        return new_hash_table
    else:
        return hash_table

''' In the first layer check for another level higher of the unobserved subclones for more intersected mutations. '''
# Input is Tree, Nodes in t=-1, Nodes in t=0, unobserved subclones before t=0
# Updates list of unobserved subclones if any new node found
def checkMoreUnobservedSubclone(Tree, prev_nodes, t0_node_list, usc_nodes_list):
    # Contains unobserved subclones between t=-1 and t=0, and nodes in t=0
    node_mut = {}
    t0_node_list.extend(usc_nodes_list)
    for n in t0_node_list:
        node_mut[n] = Tree[n].mutations

    hash_table = build_hash_table(node_mut) # Hash table should be from all nodes in t1
    updated_hash_table = eliminate_usc_sets(Tree, t0_node_list, hash_table)

    # Sort hash table to have the highest children as the first item
    sorted_hash_table = dict(sorted(updated_hash_table.items(), key=lambda items:len(items[1]), reverse=True))
    recheck_hash_table = checkHashTable(sorted_hash_table)
    usc_nodes = select_unobservedSubclones(Tree, [0], recheck_hash_table, 0, 1, usc_nodes_list)
    print("Unobserved subclones in t=0 after 1st iteration ",usc_nodes)
    #usc_nodes_list.extend(usc_nodes)
    #prev_nodes.extend(usc_nodes)
    return usc_nodes

''' List unobserved subclones by building the Hash table. '''
# Input is Tree, nodes (not unobserved) in each timepoint, already connected nodes with same set of mutations
def find_unobservedSubclone(Tree, tp_nodes, connectedNodes):
    tp_usc_nodes = {} # Unobserved subclones in each time point
    print(" Before finding unobserved subclones ",tp_nodes)
    for t1, n1 in tp_nodes.items():
        if t1+1 not in tp_nodes:
            continue
        # Update this list with unobserved nodes and t1 nodes
        t1_nodes_list = copy.deepcopy(n1)
        print(" Timepoint ",str(t1))
        # Timepoint for unobserved subclones. This will only be used to check parallel mutations later.
        usc_timepoint = t1+1

        t2_nodes = copy.deepcopy(tp_nodes[t1+1])
        node_mut = nodeMutations(Tree, t2_nodes, connectedNodes)
        hash_table = build_hash_table(node_mut) # hash table is from nodes in next timepoint
        iterNo = 1
        # After eliminating unobserved subclones updated_hash_table can be {}. Then stop here
        while(t1_nodes_list != []):
            print("Iter ",iterNo)
            print(" Nodes in t1 ",t1_nodes_list)
            updated_hash_table = eliminate_usc_sets(Tree, t1_nodes_list, hash_table)
            print(" Updated Hash table ",updated_hash_table)
            if updated_hash_table == {}: # 1st step where no more unobserved subclone is left
                break
            # Sort hash table to have the highest children as the first item
            sorted_hash_table = dict(sorted(updated_hash_table.items(), key=lambda items:len(items[1]), reverse=True))
            print(" Sorted hash table ",sorted_hash_table)
            recheck_hash_table = checkHashTable(sorted_hash_table)
            print(" After rechecking hash table ",recheck_hash_table)
            # Get the updated hash table for next iteration
            usc_nodes = select_unobservedSubclones(Tree, t1_nodes_list, recheck_hash_table, usc_timepoint, iterNo, [])
            print("usc_nodes ",usc_nodes)
            if usc_nodes == []: # 2nd step where no more unobserved subclone is left
                break

            t1_nodes_list.extend(usc_nodes)
            if iterNo == 1 and t1 == -1 and usc_nodes != []:
                new_usc_nodes = checkMoreUnobservedSubclone(Tree, t1_nodes_list, t2_nodes, usc_nodes)
                usc_nodes.extend(new_usc_nodes)
                t1_nodes_list.extend(new_usc_nodes)

            if t1+1 in tp_usc_nodes:
                tp_usc_nodes[t1+1].extend(usc_nodes)
            else:
                tp_usc_nodes[t1+1] = usc_nodes

            iterNo = iterNo+1
    print("Inside find_unobservedSubclone ",tp_nodes)
    return Tree, tp_usc_nodes

''' Look for node v such that the mutations of v that are not
    existent in u is the least, whereas v is a node in
    time point k, or the unobserved nodes between time points k and k + 1.'''
# Input is Tree, nodes in previous timepoint and mutations of not connected node in current timepoint
def select_nodeV_withLeastMutationLoss(Tree, t1_nodes, mut):
    v_mut_diff = {}
    print(mut)
    for n in t1_nodes:
        if Tree[n].pID == -2:
            continue
        #print(n," ",Tree[n].mutations)
        diff_mut = set(Tree[n].mutations) - set(mut)
        #diff_mut = set(mut) - set(Tree[n].mutations)
        v_mut_diff[n] = len(diff_mut)

    # Select the node with minimum back mutations
    print("Dict to select node v ",v_mut_diff)
    nodeV = min(v_mut_diff,key=v_mut_diff.get)
    print("Selected node v ",nodeV)
    return nodeV

''' Select the node with maximum intersections of mutations. '''
# Input is Tree, nodes in previous timepoint and mutations of the node we want to connect
def selectNode_maxCommonMutation(Tree, t1_nodes, mut):
    node_mut = {}
    for n in t1_nodes:
        if Tree[n].pID == -2:
            continue
        if set(Tree[n].mutations).issubset(set(mut)):
            node_mut[n] = len(set(Tree[n].mutations) & set(mut))

    if node_mut == {}:
        return -1
    else:
        selected_node = max(node_mut,key=node_mut.get)
        print("Selected node to connect ",selected_node)
        return selected_node

''' Finally we connect each remaining node in time $k+1$ by looking for a parent node in time $k$ 
whose consensus genotype is a subset of the node in time $k+1$. 
If such a parent node does not exist, we connect the remaining node with a node in time $k$ or 
the identified unobserved nodes in between time $k$ and $k+1$ such that the resulting back mutations is the least. '''
# Input is Tree, nodes in each timepoint, unobserved subclones between timepoints
def connect_remainingNodes(Tree, tp_nodes, tp_usc_nodes):
    for i in range(1,len(Tree)):
        if i == 7:
            print(" For Node 7 ")
        #if Tree[i].type == "unobserved":
        #    continue
        # check for nodes that are not connected
        if Tree[i].pID == -1: 
            # Check only nodes in previous timepoint and not their unobserved subclones
            print(" tp_usc_nodes ",tp_usc_nodes," tp_nodes ",tp_nodes)
            prev_tp_nodes = tp_nodes[Tree[i].timepoint-1]
            print(" List of prev_tp_nodes ",prev_tp_nodes)
            # First check for subset of mutations
            parentNode = selectNode_maxCommonMutation(Tree, prev_tp_nodes, Tree[i].mutations)
            if parentNode != -1:
                Tree[i].pID = parentNode
                Tree[Tree[i].pID].children.append(i)
            else: # Check 2nd condition
                # Unobserved nodes in this timepoint has the same timepoint
                if Tree[i].timepoint in tp_usc_nodes: # Add unobserved subclones to check
                    prev_tp_nodes.extend(tp_usc_nodes[Tree[i].timepoint])

                print(" List of prev tp nodes ",prev_tp_nodes)
                print(" Selecting node v for node ",i)
                Tree[i].pID = select_nodeV_withLeastMutationLoss(Tree, prev_tp_nodes, Tree[i].mutations)
                Tree[Tree[i].pID].children.append(i)
    return Tree

# Input child node id and parent node id
# Output is an incoming edge of type string
def getEdges(Tree):
    for i in range(len(Tree)):
        Tree[i].edges = str(Tree[i].pID)+"_"+str(i)
        #p = Tree[i].pID
        #if i not in Tree[p].children:
        #    Tree[p].children.append(i)
    return Tree
                
#tp_cluster_cells = {0: {0: [0, 1, 2], 1: [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], 2: [22], 3: [23], 4: [24, 25, 26]}, 1: {0: [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40], 1: [41, 42, 44, 45, 46, 47, 48], 2: [43], 3: [49, 50, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62], 4: [51, 52]}, 2: {0: [63], 1: [64, 66, 67, 68, 69, 70, 71, 72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89], 2: [65, 74, 75]}}
#tp_cluster_genotype = {0: {0: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1: [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0], 2: [1, 3, 1, 0, 0, 1, 1, 0, 0, 0, 0, 3, 3, 0, 3, 1, 0, 1, 3, 0], 3: [1, 1, 1, 0, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 1, 1, 0, 0], 4: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0]}, 1: {0: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1: [1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0], 2: [0, 0, 1, 0, 0, 3, 1, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 3, 0, 0], 3: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0], 4: [1, 1, 1, 1, 1, 1, 1, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0]}, 2: {0: [1, 1, 1, 1, 1, 3, 1, 0, 0, 0, 0, 1, 1, 0, 0, 3, 0, 0, 3, 0], 1: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1], 2: [1, 1, 1, 1, 0, 0, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]}}
#tp_cluster_cells = {0: {0: [0, 1, 2], 1: [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21, 22, 23], 2: [24, 25, 26]}, 1: {0: [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40], 1: [41, 42, 44, 45, 46, 47, 48, 43], 2: [49, 50, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62, 51, 52]}, 2: {0: [63, 64, 66, 67, 68, 69, 70, 71, 72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89, 65,74, 75]}}
#tp_cluster_genotype = {0: {0: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 0, 0, 0, 0, 0, 0], 1: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0], 2: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0]}, 1: {0: [1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0], 1: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0], 2: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 0]}, 2: {0: [1, 1, 1, 1, 1, 1, 1, 0, 0, 0, 0, 1, 1, 1, 1, 1, 1, 1, 1, 1]}}

#Tree, clone_node, tp_nodes = buildTree(tp_cluster_cells, tp_cluster_genotype)
#printTree(Tree)
#print(tp_nodes)
#Tree, connectedNodes = connectSameNodes(Tree)
#printTree(Tree)
#Tree, tp_usc_nodes = find_unobservedSubclone(Tree, tp_nodes, connectedNodes)
#print("Unobserved subclones in each timepoint ",tp_usc_nodes)
#Tree, connect_remainingNodes(Tree, tp_nodes, tp_usc_nodes)
#printTree(Tree)

