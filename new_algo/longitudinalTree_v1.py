import time
import argparse
from collections import OrderedDict
import itertools
from itertools import combinations, chain
import sortedcontainers
from sortedcontainers import SortedList, SortedSet, SortedDict
import random
import gc 
import copy

''' Arrange the subclones in a graph using this class. '''
class graphnode():
    def __init__(self, id_):
        # id here is the subclone id
        self.id = id_
        self.mutations=[]
        # the edge is the one above the node; edge ID is the same as node id.
        self.edge_length = []
        # root's parent ID is -1
        self.parentID = []
        self.timepoint_ = -1
        self.type = 'subclone' # can be subclone or unobserved subclone.
        #self.perc = -1
        #self.if_leaf = 1 # 1 indicates it is a leaf, 0 indicates not a leaf and -1 indicates dead.
    def getID(self):
        return self.id
    def getTimepoint(self):
        return self.timepoint_
    #def getPerc(self):
    #    return self.perc
    def getEdgeLength(self):
        return self.edge_length

class mut_Graph:
    # Constructor
    def __init__(self):
        # self.num_of_nodes = []
        self.num_of_nodes = set() # No. of mutations.
        # list representations of a graph
        self.list_of_edges = []

    def num_nodes(self):
        return self.num_of_nodes
        #return len(self.num_of_nodes)

    # Add edge to a graph
    def add_edge(self, node1, node2):
        # Add the edge from node1 to node2
        self.list_of_edges.append([node1, node2])
        self.num_of_nodes.add(node1)
        self.num_of_nodes.add(node2)

    # Print a graph representation
    def print_edge_list(self):
        print("edges: ", self.list_of_edges)
        num_of_edges = len(self.list_of_edges)
        for i in range(num_of_edges):
            print("edge ", i+1, ": ", self.list_of_edges[i])

''' Arrange the subclones in a graph using this function. '''
def inputToGraph(Gprime_file):
    file = open(Gprime_file,"r")
    line = file.readline().rstrip('\n')
    Graph = []
    Graph.append(graphnode(0)) # Append the root node.
    Graph[0].id = 0
    Graph[0].timepoint_ = 1
    Graph[0].parentID = [-1]
    node = 1
    node_tp = SortedDict()
    node_tp[1] = [0]
    while(line != ""):
        line_a = line.split('\t')
        Graph.append(graphnode(node))
        Graph[node].id = node
        Graph[node].timepoint_ = (line_a[0].split("_"))[1]
        if int(Graph[node].timepoint_) in node_tp:
            node_list = node_tp[int(Graph[node].timepoint_)]
            node_list.append(node)
            node_tp[int(Graph[node].timepoint_)] = node_list
        else:
            node_list = []
            node_list.append(node)
            node_tp[int(Graph[node].timepoint_)] = node_list

        mut_arr = []
        for i in range(1,len(line_a)):
            if int(line_a[i]) == 1:
                mut_arr.append(i-1)
        Graph[node].mutations = mut_arr
        node = node+1
        line = file.readline().rstrip('\n')
    
    for j in range(len(Graph)):
        print(" ID ",Graph[j].id," Timepoint ",Graph[j].timepoint_," Mutations ",Graph[j].mutations," Type ",Graph[j].type)
    print("Timepoints ",node_tp)
    return Graph, node_tp

''' Adds filtered unobserved subclones. Updates edges between usc and nodes in all timepoints. '''
def enroll_usc(Graph, node_s, s_pIDs):
    new_nodeId = len(Graph)
    for node, usc in node_s.items():
        usc = list(int(u) for u in usc)
        # Check here if node is single and also check mutations are same with timepoints.
        nodeIds = node.split(',')
        print(" Enroll USC nodes ",nodeIds," len ",len(nodeIds))
        if len(nodeIds) == 1:
            nid = int(nodeIds[0])
            nid_mut = Graph[nid].mutations
            print(nid_mut," ============ ",usc)
            if set(nid_mut) == set(usc):
                continue

        #print(" USC ",usc)
        Graph.append(graphnode(new_nodeId))
        Graph[new_nodeId].id = new_nodeId
        Graph[new_nodeId].mutations = usc
        Graph[new_nodeId].type = 'unobserved subclone'

        usc_key = ','.join(map(str,usc))
        usc_pid = int(s_pIDs[usc_key])
        #print(" Usc pID ",usc_pid)
        pID = []
        pID.append(usc_pid) # Get the parent node Id
        Graph[new_nodeId].parentID = pID
        p_mut = Graph[usc_pid].mutations
        edge_len = []
        edge_len.append(len(set(usc) - set(p_mut)))
        Graph[new_nodeId].edge_length = edge_len # New mutations in the usc when compared with node in t.

        for nid in nodeIds: # Connect nodes in t+1 with the unobserved subclones
            nid = int(nid)
            t2 = int(Graph[nid].timepoint_)
            #print(" Time t2 ",t2)
            node2_pIds = Graph[nid].parentID
            node2_pIds.append(new_nodeId)
            Graph[nid].parentID = node2_pIds # Append new parentIDs.

            node2_mut = Graph[nid].mutations
            node2_edge = Graph[nid].edge_length
            node2_edge.append(len(set(node2_mut) - set(usc)))
            Graph[nid].edge_length = node2_edge # New mutations in node t+1 compared to usc.

        Graph[new_nodeId].timepoint_ = float((t2 + (t2 - 1))/2) # Update usc timepoints.
        new_nodeId = new_nodeId+1

    return Graph

''' Find the largest nodes list that are connected to the unobserved subclones. '''
def find_mostConnectedNodes(nodes_list):
    node_sets={frozenset(nodes) for nodes in nodes_list}  
    updated_Nodes_list = []
    for nodes in node_sets:
        if any(nodes < ns for ns in node_sets):
            continue
        else:
            updated_Nodes_list.append(list(nodes))  
    #print(" Updated node list ",updated_Nodes_list)
    return updated_Nodes_list

''' Filter the unobserved subclones by finding node s and set A_bar. '''
def find_nodeS(usc_children):
    node_s = {} # The list of keys can be treated as nodes in set A that are connected to their corresponding usc which is value here.
    s_pIDs = {}
    maxNodes = find_mostConnectedNodes(list(usc_children.values())) # Get the maximum children or node list.

    # Find node s by selecting unobserved subclone with max no. of children. In case of tie select one with max no. of mutations.
    for subset, node in usc_children.items():
        for mnode in maxNodes:
            #print(set(node)," ",set(mnode))
            if not set(node) == set(mnode):
                continue

            parent_node = (subset.split("_"))[1]
            subsets = ((subset.split("_"))[0]).split(',')
            subsets_mut = len(subsets)
            #print(subset," ",node_count," ",subsets_mut," ",node)
            key_node = ','.join(map(str,node))
            s_key = ','.join(map(str,subsets))
            if key_node in node_s:
                mut_count = len(node_s[key_node]) 
                if subsets_mut > mut_count: # In case of tie choose the node with max no. of mutations
                    node_s[key_node] = subsets
                    s_pIDs[s_key] = parent_node
            else:
                node_s[key_node] = subsets
                s_pIDs[s_key] = parent_node

    print(" Node s ",node_s) 
    return node_s,s_pIDs

''' Find all subsets first and then get proper subsets.''' 
def find_subsets(set1, set2): # set1 is mutations at time t. set2 is mutations at time t+1.
    proper_subsets = set()
    iterList = list(set2)
    # We started at len = set1 + 1 and went till len of set2 - 1 such that proper subsets and proper supersets are obtained. 
    allsubsets = chain.from_iterable(combinations(iterList, ss) for ss in range(len(set1)+1, len(iterList)))
    #print(" allsubsets ",allsubsets)

    for als in allsubsets: # From all subsets only choose the proper subsets.
        #print(" als ",als)
        if set(set1) < set(als):
            proper_subsets.add(als)

    gc.collect()
    return proper_subsets

''' Find the common subset of mutations(omega set) shared among all remaining nodes in t+1. '''
def find_commonSubsets(all_n2_mutations):
    print(" All n2 mutations ",all_n2_mutations)
    omega_set = []
    for i in range(0,len(all_n2_mutations)):
        for j in range(1,len(all_n2_mutations)):
            if i == j:
                continue
            if len(omega_set) > 0:
                for k in range(0, len(omega_set)):
                    temp_set = set(omega_set[k]) & set(all_n2_mutations[j])
                    if not len(omega_set) == 0:
                        omega_set.append(temp_set)
                        omega_set.remove(omega_set[k])
            else:
                temp_set = set(all_n2_mutations[i]) & set(all_n2_mutations[j])
                if not len(temp_set) == 0:
                    omega_set.append(temp_set)
    print(" Omega set ",omega_set)
    return omega_set

''' Find the proper subset of the nodes in time t that contain all mutations in \omega. '''
def find_ps_tNodes(Graph, timepoint, omega_set):
    print("====== Inside find PS t nodes ==========")
    proper_subsets = []
    for i in range(len(Graph)):
        if not int(Graph[i].timepoint_) == int(timepoint):
            continue
        #print(" Node id ",i)
        mut = Graph[i].mutations
        print(" Time t mutations ",mut)
        if len(mut) > 0 and len(omega_set) > 0 and set(mut) < omega_set[0]:
            print(" Proper subset ",mut)
            proper_subsets.append(find_subsets(mut, omega_set[0]))
        if len(mut) == 0 and len(omega_set) > 0:
            proper_subsets.append(omega_set[0])
    print(" Proper subset after checking ",proper_subsets)
    return proper_subsets

''' Check if the unobserved subclones is a superset of mutations in t. '''
def find_supersets_tNodes(Graph, timepoint, potential_usc_sets):
    print("====== Inside find superset of t nodes ==========")
    t_supersets = []
    for i in range(len(Graph)):
        if not int(Graph[i].timepoint_) == int(timepoint):
            continue
        #print(" Node id ",i)
        t_mut = Graph[i].mutations
        #print(" Time t mutations ",t_mut)
        for mut in potential_usc_sets:
            #print(" Mutations in USC ",mut)
            if set(mut) > set(t_mut):
                if mut not in t_supersets:
                    t_supersets.append(mut)
    #print(" Supersets ",t_supersets)
    return t_supersets

''' For each cluster in the assignment file get the cells. cluster --> [cellIDs]'''
def get_cell_cluster(assignment):
    with open(assignment,'r') as f:
        lines = f.readlines()

    assignment_list = ((lines[1].split('\t'))[2]).split(" ")
    cluster_cells = {} # key is clusterID, value is cell IDs.
    #print(len(assignment_list))
    cell_count = 0
    for i in assignment_list:
        if '\n' in i:
            i = i.replace("\n","")
        if i in cluster_cells:
            temp_cell_list = cluster_cells[i]
            temp_cell_list.append(cell_count)
            cluster_cells[i] = temp_cell_list
            cell_count = cell_count+1
        else:
            temp_cell_list = []
            temp_cell_list.append(cell_count)
            cluster_cells[i] = temp_cell_list
            cell_count = cell_count+1
    return cluster_cells

''' From the input Gprime matrix get the mutations for each cluster. '''
def get_cluster_mutations(Gprime_file):
    file = open(Gprime_file,"r")
    line = file.readline().rstrip('\n')
    cluster_mutations = {}
    while(line != ""):
        line_a = line.split('\t')
        clusterID = (line_a[0].split("_"))[0]
        if clusterID in cluster_mutations:
            mut_arr = cluster_mutations[clusterID]
        else:
            mut_arr = []
        for i in range(1,len(line_a)):
            if int(line_a[i]) == 1:
                if (i-1) not in mut_arr:
                    mut_arr.append(i-1)

        cluster_mutations[clusterID] = mut_arr
        line = file.readline().rstrip('\n')
    return cluster_mutations

''' Get a dictionary where mutation is the key and cells are the values. '''
def get_mutation_cells(cluster_cells, cluster_mutations, all_n2_mut):
    mutation_cells = {}
    n2_mut = []
    for anm in all_n2_mut:
        for sub_anm in anm:
            if sub_anm not in n2_mut:
                n2_mut.append(sub_anm)
    #print(" All n2 mut ",n2_mut)
    for cluster, cells in cluster_cells.items():
        mut_arr = cluster_mutations[cluster]
        for mut in mut_arr:
            if mut not in n2_mut:
                continue
            if mut in mutation_cells:
                cell_arr = mutation_cells[mut]
                cell_arr.extend(cells)
                mutation_cells[mut] = cell_arr
            else:
                mutation_cells[mut] = cells
    return mutation_cells

''' Find the connected components using the union-find algorithm. '''
def find(parent, x):
    if x != parent[x]:
        parent[x] = find(parent, parent[x])
    return parent[x]

def union(parent, x, y):
    parent_x = find(parent, x)
    parent_y = find(parent, y)
    if parent_x != parent_y:
        # setting parent of x as y
        parent[parent_y] = parent_x

''' The parent list where mutations are grouped by similar edges. Process that list to find possible unobserved subclones. '''
def find_possible_usc(parent,parent_dict):
    usc_nodes = {}
    set_parent = set(parent)
    for sp in set_parent:
        temp_list = []
        for i in range(len(parent)):
            if parent[i] == sp:
                temp_list.append(i)
            usc_nodes[sp] = temp_list

    for k,v in usc_nodes.items(): # Update the list with the mutation values instead of indices.
        mut_list = v
        for mut, index in parent_dict.items():
            for j in range(len(mut_list)):
                if mut_list[j] == index:
                    mut_list[j] = mut
        usc_nodes[k] = mut_list
    return usc_nodes

''' Use this parent dict to use the union-find algorithm. '''
def getParent_dict(parent):
    parent_dict = {}
    for i in range(len(parent)): # Create an index for each parent where node is key and index is value.
        parent_dict[parent[i]] = i
    return parent_dict

''' If any unconnected mutations still have number of cells greater than 2 then include them as single usc. '''
def unconnected_mutations(parent,mutation_cells,all_n2_mutations):
    #mutation_cells = {1: [1,2], 2: [3,4,5], 3: [3]}
    #parent = [4,5]
    mut_set = set(mutation_cells.keys())
    unconnected_mut = set()
    mut_usc = set()
    #print(" Mutation set ",len(mut_set))
    #print(" Parent ",len(parent))
    for m in mut_set:
        if m not in parent:
            um_in_all_nodes = 0
            for anm in all_n2_mutations: # The unconnected mutation should also appear in all subclones in t+1.
                if m in anm:
                    um_in_all_nodes = um_in_all_nodes+1
            if um_in_all_nodes == len(all_n2_mutations): 
                unconnected_mut.add(m)
    for um in unconnected_mut:
        mut_cells = mutation_cells[um]
        if len(mut_cells) >= 2:
            mut_usc.add(um)

    print("Unconnected Mutation usc ",mut_usc)
    #print(" Unconnected mutations ",unconnected_mut)
    return mut_usc

''' Build the mutation graph where each node represents a mutation and each node is connected if they have same list of cells. 
After grouping mutations using union-find algorithm return the possible unobserved subclones. '''
def process_mut_Graph(mutation_cells,all_n2_mutations):
    m_graph = mut_Graph()
    for mut1, cells1 in mutation_cells.items(): # Build the mutation graph.
       for mut2, cells2 in mutation_cells.items():
           if mut1 == mut2 or len(cells1) != len(cells2):
               continue
           cells1.sort()
           cells2.sort()
           if cells1 == cells2:
               m_graph.add_edge(mut1, mut2) 
    # Mutation graph built with connected edges.
    #m_graph.print_edge_list()
    parent = list(m_graph.num_nodes())
    parent_dict = getParent_dict(parent) # Create a dict to map mutations with the index array.
    #print(parent_dict)
    #print("Actual parent: ", parent) # Modify the union-find algorithm to keep the mutations. Think how to do that.
    parent_list = list(parent_dict.values())
    #print("Before processed parent: ",parent_list)
    # for each edge a,b check if a is connected to b or not
    for a,b in m_graph.list_of_edges: # Use union-find algorithm to group mutations 
        #print(a," ",b)
        a = parent_dict[a] # Pass the index values for union-find algorithm.
        b = parent_dict[b]
        union(parent_list, a, b)

    #print("Processed parent: ", parent_list) # This parent list now has mutations grouped which is replaced by its parent.
    # Reverse parent_list by replacing the values with actual mutations.
    for mut, index in parent_dict.items():
        for i in range(len(parent_list)):
            if parent_list[i] == index:
                parent_list[i] = mut
    #print("Replaced parent: ",parent_list)

    unconn_mut = unconnected_mutations(parent,mutation_cells,all_n2_mutations)
    usc_nodes_dict = find_possible_usc(parent_list,parent_dict)
    usc_nodes_list = list(usc_nodes_dict.values())
    if unconn_mut != {}: # Add the single unconnected mutations as usc.
        for um in unconn_mut:
            temp_list = []
            temp_list.append(um)
            usc_nodes_list.append(temp_list)

    print(" USC nodes ",usc_nodes_list)
    return usc_nodes_list

''' Filter the mutations that are already obtained from step 2 before proceeding to step 3. '''
def filter_mutations(ps_fromT, mutation_cell):
    mutation_set = set(mutation_cell.keys())
    #print(" Mutation set ",mutation_set)
    existing_mut = set()
    for ps in ps_fromT:
        for sub_ps in ps:
            #print(" sub_ps ",sub_ps)
            if sub_ps in mutation_set:
                del mutation_cell[sub_ps]
                existing_mut.add(sub_ps)
    #print(" Existing mutation ",existing_mut)
    #print(" Mutation cell ",mutation_cell.keys())
    return mutation_cell, existing_mut

''' Get the usc_children dict where usc and its parent node is the key and nodes in time t+1 are its children. '''
def get_usc_children(Graph, timepoint, possible_usc_list, usc_children):
    usc_keys = set()
    for i in range(len(Graph)): # First check if it is a superset of nodes in time t.
        if not int(Graph[i].timepoint_) == int(timepoint):
            continue
        t_mut = Graph[i].mutations
        for pusc in possible_usc_list:
            if set(t_mut).issubset(set(pusc)):
                key = ",".join(map(str,pusc))+"_"+str(i)
                usc_keys.add(key)

    for j in range(len(Graph)): # Then check if it is a subset of nodes in time t+1 and connect them.
        if not int(Graph[j].timepoint_) == int(timepoint+1):
            continue
        t1_mut = Graph[j].mutations
        for pusc in possible_usc_list:
            #print(" PUSC ",pusc," t1 mut ",t1_mut)
            if set(pusc).issubset(set(t1_mut)):
                for uk in usc_keys:
                    uk_usc = (uk.split('_')[0]).split(',')
                    uk_usc = (int(i) for i in uk_usc) # Convert the mutations to int from str.
                    #print(" USC key ",uk_usc)
                    if set(uk_usc) == set(pusc):
                        if uk in usc_children:
                            temp_list = usc_children[uk]
                            temp_list.append(j)
                            usc_children[uk] = temp_list
                        else:
                            temp_list = []
                            temp_list.append(j)
                            usc_children[uk] = temp_list

    print(" Time ",timepoint," USC keys ",usc_keys)     
    print(" USC children ",usc_children)
    return usc_children

''' Build the hash table with nodes and mutations. '''
def build_hash_table(all_n2_mutations):
    mut_dict = {}
    hash_table = {}
    for i in range(len(all_n2_mutations)):
        mut_list = all_n2_mutations[i]
        for mut in mut_list:
            if mut in mut_dict:
                temp_list = mut_dict[mut]
                temp_list.append(i)
                mut_dict[mut] = temp_list
            else:
                temp_list = []
                temp_list.append(i)
                mut_dict[mut] = temp_list

    for mut, nodes in mut_dict.items():
        hash_table_key = ",".join(map(str,nodes))
        if hash_table_key in hash_table:
            temp_list = hash_table[hash_table_key]
            temp_list.append(mut)
            hash_table[hash_table_key] = temp_list
        else:
            temp_list = []
            temp_list.append(mut)
            hash_table[hash_table_key] = temp_list
    print(" Hash table ",hash_table)
    return hash_table
 
''' Second, before you enumerate all the proper subset of the nodes in time t, first find the common subset of mutations shared among all remaining nodes in t+1 (defined as X above). Let \omega this common subset of mutations. Then find the proper subset of the nodes in time t that contain all mutations in \omega. '''
''' Finding unobserved subclones by looking for proper subsets and supersets between two timepoints. '''
def find_unobservedSubclones(Graph, node_tp, cluster_cells, cluster_mutations):
    usc_children = {}
    for time1, nodes in node_tp.items():
        t1_nodes = nodes
        if time1+1 in node_tp: # Compare nodes in t with nodes in t+1.
            t2_nodes = node_tp[time1+1]
        else:
            continue

        print(" Timepoint ==== ",time1)
        #usc_timepoint = float((time1+(time1+1))/2) # Timepoint for the unobserved subclones.
        #print(" Usc timepoint ",usc_timepoint)
        connected_n2 = []
        all_n2_mutations = []

        for n1 in t1_nodes:
            n1_mut = Graph[n1].mutations
            #print(" Node 1 ",n1," N1 mut ",n1_mut)
            start = time.time()
            
            for n2 in t2_nodes: # Compare mutations of nodes between two timepoints to find unobserved subclones.
                if n2 in connected_n2:
                    continue
                n2_mut = Graph[n2].mutations
                if n1_mut == n2_mut:
                    pID_list = []
                    pID_list.append(n1)
                    Graph[n2].parentID = pID_list
                    Graph[n2].edge_length = [0] # since there are no new mutations, the number is 0.

                    connected_n2.append(n2) # Mark connected nodes in t and t+1.
                    if n2_mut in all_n2_mutations: # Second check to see if same nodes exist in t and t+1 then don't find their supersets.
                        all_n2_mutations.remove(n2_mut)
                    continue

                if n2_mut not in all_n2_mutations:
                    all_n2_mutations.append(n2_mut)
            print(" Length of Graph ",len(Graph))

        # Build the hash table here
        hash_table = build_hash_table(all_n2_mutations)

        #mutation_cell = get_mutation_cells(cluster_cells, cluster_mutations, all_n2_mutations)
        #print(" Mutation cells ",mutation_cell.keys())

        #omega_set = find_commonSubsets(all_n2_mutations) # This is step 2
        #ps_tNodes = find_ps_tNodes(Graph, time1, omega_set) # This is also step 2
        #ps_tNodes = [[100]]
        #print(" Step 2 Proper subset ",ps_tNodes)
        #if ps_tNodes != []:
            # Before filtering check if existing_mut has no. of cells >= 2 in  mutation_cell. If yes, then keep them.
            #keep_existing_mut = []
            #for mut, cell in mutation_cell.items():
            #    for ps in ps_tNodes:
            #        for sub_ps in ps:
            #            if int(sub_ps) == int(mut):
            #                if len(cell) >= 2:
            #                    if sub_ps not in keep_existing_mut:
            #                        keep_existing_mut.append(sub_ps)

            #print(" Keep existing mutations ",keep_existing_mut)
            #mutation_cell, existing_mut = filter_mutations(ps_tNodes, mutation_cell)
            #potential_usc_sets = process_mut_Graph(mutation_cell,all_n2_mutations) ## Step 3 starts here
            #updated_n2_mutations, existing_mut = filter_cs_mutations(ps_tNodes,all_n2_mutations)
            #print(" Proper subset not empty so existing mut ",existing_mut)
            
            #hash_table = build_hash_table(updated_n2_mutations)
            #potential_usc_sets = list(hash_table.values())

            #if len(existing_mut) > 0:
                #print(" Here in the condn ")
                #print(" Potential USC ",potential_usc_sets)
            #    for i in range(len(potential_usc_sets)):
            #        potential_usc_sets[i].extend(list(existing_mut))
            #if keep_existing_mut != []:
            #    potential_usc_sets.append(keep_existing_mut)
            #print(" Potential usc sets ",potential_usc_sets)
            #print(" New logic ends here ... ")
       # else:
            #potential_usc_sets = process_mut_Graph(mutation_cell,all_n2_mutations)
        potential_usc_sets = list(hash_table.values())
        if time1 == 1:
            final_possible_usc = potential_usc_sets # At t=1 or when root we directly use the sets we get without looking for supersets. 
            print(" Final possible usc ",final_possible_usc)
        else:
            final_possible_usc = find_supersets_tNodes(Graph, time1, potential_usc_sets) 
            print(" Final possible usc ",final_possible_usc)

        usc_children = get_usc_children(Graph, time1, final_possible_usc, usc_children) # Update the usc node.
        print(" USC children ",usc_children)

    #print(" Timepoint usc dict ",len(timepoint_usc[2]))
    return Graph, usc_children

''' For such nodes u we find a node v in time k
that has a set of mutations that are also in u, and such a set’s size is the largest.
If there are more than one node v then we select one randomly. '''
def find_nodeV(Graph, node_tp):
    timepoints = node_tp.keys()
    for i in range(1,len(timepoints)):
        t2 = timepoints[i]
        t1 = t2 - 1
        print(" Time ",timepoints[i])
        t2_nodes = node_tp[t2]
        t1_nodes = node_tp[t1]

        for u in t2_nodes:
            if not Graph[u].parentID == []:
                continue
            u_mut = Graph[u].mutations
            v_mutations = []
            v_nodes = []
            for v in t1_nodes:
                v_mut = Graph[v].mutations
                if set(v_mut).issubset(set(u_mut)):
                    v_mutations.append(v_mut)
                    v_nodes.append(v)

            largest_v_mut = find_mostConnectedNodes(v_mutations) # Find the largest mutation set size.
            print(" largest_v_mut ",largest_v_mut)
            if len(largest_v_mut) > 1:
                random_index = random.randrange(len(largest_v_mut)) # Find a random mutation set if more than one exist.
                largest_v_mut = largest_v_mut[random_index]
            else:
                largest_v_mut = largest_v_mut[0]

            for i in range(len(v_mutations)):
                if set(largest_v_mut) == set(v_mutations[i]):
                    v_node_index = i
                    break

            pID = []
            pID.append(v_nodes[v_node_index])
            Graph[u].parentID = pID

            edge_len = []
            edge_len.append(len(set(u_mut) - set(largest_v_mut)))
            Graph[u].edge_length = edge_len

    return Graph

''' Save the graph. '''
def save_graph(G, out_f):
    f = open(out_f, "w")
    #f.write("\t".join(['ID','Timepoint','PID','Type','Mutations','Edge length'])+ "\n")
    for i in range(len(G)):
        #print(",".join(map(str,G[i].parentID)))
        f.write("\t".join([str(G[i].id), str(G[i].timepoint_), ",".join(map(str,G[i].parentID)), str(G[i].type), ",".join(map(str,G[i].mutations)), ",".join(map(str,(G[i].edge_length)))]) + "\n")
    f.close()

''' Get timepoints for each cell from the simulation file. '''
def get_cell_timepoints(celltp_file):
    file = open(celltp_file,"r")
    line = file.readline().rstrip('\n')
    cell_tp = {}
    while(line != ""):
        line_a = line.split('\t')
        #print(line_a)
        tp = (line_a[0].split('_'))[0]
        cells = line_a[1].split(';')
        for cell in cells:
            cell_tp[cell] = tp
        line = file.readline().rstrip('\n')
    return cell_tp

''' Get # of cells in each cluster. '''
def getCellsOfEachCluster(cell_tp, cluster_cell):
    cc_time = {}
    cc_len_time = {}
    for cluster,cells in cluster_cell.items():
        for cell in cells:
            celltp = cell_tp[str(cell)]
            key = cluster+'_'+celltp
            if key in cc_time:
                temp_list = cc_time[key]
                temp_list.append(cell)
                cc_time[key] = temp_list
            else:
                temp_list = []
                temp_list.append(cell)
                cc_time[key] = temp_list
    for cluster,cells in cc_time.items():
        cc_len_time[cluster] = len(cells)
    print(" Cells in each cluster ",cc_len_time)

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (Input to the Longitudinal tree)")
parser.add_argument("-cc", "--cc",dest ="cc", help="Assignment file indicating cells belonging to clusters.")
parser.add_argument("-celltp", "--celltp",dest = "celltp", help="Cell timepoints.")
args = parser.parse_args()
Graph, node_tp = inputToGraph(args.input) # Creates a graph from the input G' matrix.

cluster_cell = get_cell_cluster(args.cc)
#print(" Cluster cell ")
#print(cluster_cell)
#print(" ============ ")
cluster_mutation = get_cluster_mutations(args.input)
print(cluster_mutation)
#print(" ============ ")
cell_tp = get_cell_timepoints(args.celltp)
getCellsOfEachCluster(cell_tp, cluster_cell)
#print(" Cell timepoints ",cell_tp)

up_Graph, usc_children = find_unobservedSubclones(Graph, node_tp, cluster_cell, cluster_mutation) # Updates the graph by including unobserved subclones, adding edges and branch length.
node_s,s_pIDs = find_nodeS(usc_children)
up_Graph = enroll_usc(Graph, node_s, s_pIDs)
print(" Node TP ",node_tp)
up_Graph = find_nodeV(Graph, node_tp)

for i in range(len(up_Graph)):
    print(up_Graph[i].id," ",up_Graph[i].timepoint_," ",up_Graph[i].parentID," ",up_Graph[i].type," ",up_Graph[i].mutations," ",up_Graph[i].edge_length)

#save_graph(up_Graph,"graph.csv")
