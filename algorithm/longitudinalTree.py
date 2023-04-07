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
import numpy as np
import statistics

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

#class mut_Graph:
    # Constructor
#    def __init__(self):
        # self.num_of_nodes = []
#        self.num_of_nodes = set() # No. of mutations.
        # list representations of a graph
#        self.list_of_edges = []

#    def num_nodes(self):
#        return self.num_of_nodes
        #return len(self.num_of_nodes)

    # Add edge to a graph
#    def add_edge(self, node1, node2):
        # Add the edge from node1 to node2
#        self.list_of_edges.append([node1, node2])
#        self.num_of_nodes.add(node1)
#        self.num_of_nodes.add(node2)

    # Print a graph representation
#    def print_edge_list(self):
#        print("edges: ", self.list_of_edges)
#        num_of_edges = len(self.list_of_edges)
#        for i in range(num_of_edges):
#            print("edge ", i+1, ": ", self.list_of_edges[i])

''' Arrange the subclones in a graph using this function. '''
def inputToGraph(Gprime_file):
    file = open(Gprime_file,"r")
    line = file.readline().rstrip('\n')
    clone_node = {} # Map the clones to the nodes in Graph for downstream analysis.
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
        clone_node[line_a[0]] = node
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
    return Graph, node_tp, clone_node

''' Build the graph after correcting the clusters based on likelihood values. '''
def build_graph(corrected_tp_clusters):
    clone_node = {} # Map the clones to the nodes in Graph for downstream analysis.
    Graph = []
    Graph.append(graphnode(0)) # Append the root node.
    Graph[0].id = 0
    Graph[0].timepoint_ = 1
    Graph[0].parentID = [-1]
    node = 1
    node_tp = SortedDict()
    node_tp[1] = [0]
    for tp, clusters in corrected_tp_clusters.items():
        for cluster, genotype in clusters.items():
            Graph.append(graphnode(node))
            Graph[node].id = node
            Graph[node].timepoint_ = tp
            clone_node[cluster+"_"+str(tp)] = node
            if int(Graph[node].timepoint_) in node_tp:
                node_list = node_tp[int(Graph[node].timepoint_)]
                node_list.append(node)
                node_tp[int(Graph[node].timepoint_)] = node_list
            else:
                node_list = []
                node_list.append(node)
                node_tp[int(Graph[node].timepoint_)] = node_list
            
            mut_arr = []
            for i in range(len(genotype)):
                if int(genotype[i]) == 1:
                    mut_arr.append(i)
            Graph[node].mutations = mut_arr
            node = node+1
    
    # print the graph with its timepoint and mutations
    for j in range(len(Graph)):
        print(" ID ",Graph[j].id," Timepoint ",Graph[j].timepoint_," Mutations ",Graph[j].mutations," Type ",Graph[j].type)
    print("Timepoints ",node_tp)
    return Graph, node_tp, clone_node

''' Get cluster genotype by taking a majority vote for each cell. '''
def vote_cluster_genotypes(cell_gs_list):
    cell_gs_len = len(cell_gs_list)
    cg_len = len(cell_gs_list[0])
    #print(" Inside Voted cluster genotype ")
    #print(cell_gs_len," ",cg_len)

    #print(cell_gs_list[0])
    #print("============================")
    voted_cg = []

    #total0s = 0
    #total1s = 0
    for j in range(0,cg_len):
        total0s = 0
        total1s = 0
        for i in range(0,cell_gs_len):
            if cell_gs_list[i][j] == '0':
                total0s = total0s+1
            if cell_gs_list[i][j] == '1':
                total1s = total1s+1
        #print(" Total 0s ",total0s," Total 1s ",total1s)
        if total0s > total1s:
            voted_cg.append('0')
        if total1s >= total0s:
            voted_cg.append('1')
    return voted_cg

''' Get corrected cluster genotype based on true cells belonging to a clone. cell_genotype is from inferred consensus genotype of BnpC. '''
def getCorrectedClusterGenotype(cell_genotype, tp_true_clusters):
    sorted_timepoints = list(tp_true_clusters.keys())
    #print(" Sorted timepoints ",sorted_timepoints)
    corrected_tp_clusters = {}
    for timepoint, clusters in tp_true_clusters.items():
        cluster_cell_genotype = {}
        for cluster, cells in clusters.items():
            for cell in cells:
                #print(" Cell starts from ",cell)
                cell_gen = cell_genotype[int(cell)]
                #print(" Cluster ",cluster," Cell gen ",cell_gen)
                if cluster in cluster_cell_genotype:
                    t_list = cluster_cell_genotype[cluster]
                    t_list.append(cell_gen)
                    cluster_cell_genotype[cluster] = t_list
                else:
                    t_list = []
                    t_list.append(cell_gen)
                    cluster_cell_genotype[cluster] = t_list
            #print("==============================")

        #print(" Cell cluster ",cluster_cell_genotype.keys())
        corrected_cluster_genotype = {}
        for cluster, cell_gs in cluster_cell_genotype.items():
            #print(" Cluster ",cluster)
            #if len(cell_gs) == 1:
            #    cluster_genotype[cluster] = cell_gs[0]
            #    continue
            voted_cg = vote_cluster_genotypes(cell_gs)
            #print(" Cluster ",cluster," Voted CG ",voted_cg)
            #print("************************")
            corrected_cluster_genotype[cluster] = voted_cg
        #print(" Voted Cluster genotype ",corrected_cluster_genotype)
        corrected_tp_clusters[timepoint] = corrected_cluster_genotype
    return corrected_tp_clusters

''' Adds filtered unobserved subclones. Updates edges between usc and nodes in all timepoints. '''
def enroll_usc(Graph, node_s, s_pIDs):
    new_nodeId = len(Graph)
    usc_nodes = {}
    for node, usc in node_s.items():
        usc = list(int(u) for u in usc)
        # Check here if node is single and also check mutations are same with timepoints.
        nodeIds = node.split(',')
        print(" Enroll USC nodes ",nodeIds," len ",len(nodeIds))
        if len(nodeIds) == 1:
            #nid = int(nodeIds[0])
            #nid_mut = Graph[nid].mutations
            #print(nid_mut," ============ ",usc)
            #if set(nid_mut) == set(usc):
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

            # If a node is already connected then we don't reconnect it any further
            node2_pIds = Graph[nid].parentID
            if node2_pIds == []:
                node2_pIds.append(new_nodeId)
            Graph[nid].parentID = node2_pIds # Append new parentIDs.

            node2_mut = Graph[nid].mutations
            node2_edge = Graph[nid].edge_length
            node2_edge.append(len(set(node2_mut) - set(usc)))
            Graph[nid].edge_length = node2_edge # New mutations in node t+1 compared to usc.

        usc_timepoint = float((t2 + (t2 - 1))/2)
        Graph[new_nodeId].timepoint_ = usc_timepoint # Update usc timepoints.
        if usc_timepoint in usc_nodes:
            temp = usc_nodes[usc_timepoint]
            temp.append(new_nodeId)
            usc_nodes[usc_timepoint] = temp
        else:
            temp = []
            temp.append(new_nodeId)
            usc_nodes[usc_timepoint] = temp
        
        new_nodeId = new_nodeId+1

    return Graph, usc_nodes

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
        #if i == 0: # Need to change this hard-coded code
        #    i == 1
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

''' Get the cells of each cluster from the simulated file. '''
def sim_cell_cluster(clusterFile):
    with open(clusterFile,'r') as f:
        lines = f.readlines()

    cell_cluster = {}
    for line in lines:
        #print(line.split('\t'))
        clusterID = ((line.split('\t')[0]).split('_'))[1]
        cells = (line.split('\t')[1].strip()).split(';')
        cell_cluster[clusterID] = cells

    return cell_cluster

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

''' Build the hash table with nodes and mutations to reduce complexity while listing the subsets. '''
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
def find_unobservedSubclones(Graph, node_tp):
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
            #print(" Length of Graph ",len(Graph))

        # Build the hash table here
        hash_table = build_hash_table(all_n2_mutations)
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
    print(" ============== Inside finding node v ============ ")
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
                #print(u_mut," ------------- ",v_mut)
                if set(v_mut).issubset(set(u_mut)):
                    #print(" True ")
                    v_mutations.append(v_mut)
                    v_nodes.append(v)

            largest_v_mut = find_mostConnectedNodes(v_mutations) # Find the largest mutation set size.
            print(" largest_v_mut ",largest_v_mut)
            if largest_v_mut == []:
                continue

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

''' Connect nodes that doesn't have same or subsets of mutations based on largest set of mutations. '''
def connect_danglingNodes(Graph, node_tp, usc_nodes):
    int_timepoints = set(node_tp.keys())
    usc_timepoints = set(usc_nodes.keys())
    all_timepoints = sorted(int_timepoints.union(usc_timepoints))
    print(" Sorted all timepoints ",all_timepoints)
    for tp in all_timepoints:
        if tp == 1:
            continue
        if tp % 1 == 0: # Only check for nodes in integer timepoints to connect them with unobserved nodes or nodes in previous timepoint.
            print(" Timepoint ",tp)
            check_usc_nodes = []
            if (tp - 0.5) in all_timepoints:
                print(" Inside usc_node check ",tp)
                check_usc_nodes = usc_nodes[tp - 0.5]
            else:
                t1_nodes = node_tp[tp-1]
            t2_nodes = node_tp[tp]
            for u in t2_nodes:
                if not Graph[u].parentID == []:
                    continue
                print(" Node ",u," has no parent ")
                print(" USC nodes ",check_usc_nodes)
                u_mutations = Graph[u].mutations
                if check_usc_nodes != []:
                    nodes_to_consider = check_usc_nodes # Check to connect with unobserved subclone if present.
                else:
                    nodes_to_consider = t1_nodes # If unobserved subclone absent then connect with nodes in previous timepoint.
                node_mut_size = {}
                node_mutations_list = {}
                #print(" Node to consider ",nodes_to_consider)
                for node in nodes_to_consider:
                    #print(" Node to consider ",node)
                    node_mutations = Graph[node].mutations
                    #mut_size = len(set(node_mutations)-set(u_mutations)) + len(set(u_mutations)-set(node_mutations))
                    mut_size = len(set(node_mutations) - set(u_mutations))
                    node_mut_size[node] = mut_size
                    node_mutations_list[node] = node_mutations
                min_node = min(node_mut_size, key=node_mut_size.get) # Found the node with minimum mutations set.
                print(" Node with minimum mutation sets ",node_mut_size)
                pID = []
                pID.append(min_node)
                Graph[u].parentID = pID # Make that unobserved subclone parent.
                edge_len = []
                edge_len.append(len(set(u_mutations) - set(node_mutations_list[min_node])))
                #print(" Edge length ",edge_len)
                Graph[u].edge_length = edge_len
    return Graph

''' Get a dictionary from the G', G'' and D matrix to get the genotypes. '''
def getGenotype_dict(genotype, gtype):
    file = open(genotype,"r")
    line = file.readline().rstrip('\n')
    genotype_dict = {}
    if gtype == "D":
        cellCount = 0
        firstLine = True
        while(line != ""):
            #if firstLine == True:
            #    line = file.readline().rstrip('\n')
            #    firstLine = False
            #    continue
            cg = ','.join(line.split('\t'))
            genotype_dict[cellCount] = cg
            cellCount = cellCount+1
            line = file.readline().rstrip('\n')
    elif gtype == "G":
        cellCount = 0
        firstLine = True
        while(line != ""):
            if firstLine == True:
                line = file.readline().rstrip('\n')
                firstLine = False
                continue
            #cg = ','.join((line.split('\t'))[1:])
            cg = (line.split('\t'))[1:]
            #print(" CG len ",len(cg)," ",cg)
            genotype_dict[cellCount] = cg
            cellCount = cellCount+1
            line = file.readline().rstrip('\n')
    elif gtype == "cluster":
        while(line != ""):
            cluster_time = (line.split('\t'))[0]
            cg = ','.join((line.split('\t'))[1:])
            genotype_dict[cluster_time] = cg
            line = file.readline().rstrip('\n')
    return genotype_dict

#''' Find the threshold to determine that a clone does not exist in a timepoint. '''
#def find_theta(all_clones_prob):
#    prob_list = list(all_clones_prob.values())
#    prob_list.sort()
#    theta = statistics.median(prob_list)
#    return theta

#''' Calculate the probability that a subclone truly exists in the timepoint. '''
#def subclone_truly_exists(cluster_true_cells, ncells_t):
#    W_clones = [] # Wrong clones that doesn't meet threshold.
#    all_clones_prob = {} # Prob of all clones at a timepoint.
#    for cluster, cells in cluster_true_cells.items():
#        noOfcells = len(cells) # No. of supporting reads r.
#        p = noOfcells / ncells_t # Success probability.
        #print(" No of cells ",noOfcells," N cells t ",ncells_t," Prob geometric ",p)
#        prob_clones = ((1 - p)**noOfcells) * p
#        print(" Cluster ",cluster," noOfCells ",noOfcells," prob that it exists ",prob_clones)
#        all_clones_prob[cluster] = prob_clones
#    theta = find_theta(all_clones_prob)
#    for cluster, prob in all_clones_prob.items():
#        if prob >= theta:
#            W_clones.append(cluster)
#    print(" Wrong clones ",W_clones)
#    return W_clones

''' Calculate the L_ik for each cells. '''
def calculate_likelihood(cell_cluster_prob):
    cluster_true_cells = {} # Note the cells supporting the log likelihood check
    cells_not_supported = {}
    for cc, prob in cell_cluster_prob.items():
        cell = cc.split('_')[0]
        cluster = cc.split('_')[1]
        P_i_from_k = prob
        P_i_not_k = 0.0
        #print(" Cell ",cell," cluster ",cluster)
        for cc1, prob1 in cell_cluster_prob.items():
            if cc == cc1:
                continue
            cell1 = cc1.split('_')[0]
            cluster1 = cc1.split('_')[1]
            if int(cell) == int(cell1) and int(cluster) != int(cluster1):
                #print(" Cell1 ",cell1," Cluster1 ",cluster1)
                P_i_not_k = P_i_not_k + prob1
                #print(" P i not from k ",P_i_not_k)
        #print("======================")
        if P_i_not_k == 0.0 or P_i_from_k == 0.0:
            L_ik = 0
        else:
            L_ik = np.log(P_i_from_k/P_i_not_k)
        #if int(cell) == 599:
        #    print(" P i from k ",P_i_from_k," P i not from k ",P_i_not_k)
        #    print(" Log likelihood ",L_ik)
        if L_ik > 0:
            if cluster in cluster_true_cells:
                cells_list = cluster_true_cells[cluster]
                cells_list.append(cell)
                cluster_true_cells[cluster] = cells_list
            else:
                cells_list = []
                cells_list.append(cell)
                cluster_true_cells[cluster] = cells_list
        else:
            if cluster in cells_not_supported:
                cells_list = cells_not_supported[cluster]
                cells_list.append(cell)
                cells_not_supported[cluster] = cells_list
            else:
                cells_list = []
                cells_list.append(cell)
                cells_not_supported[cluster] = cells_list
            
    #print(" Clusters true cells ",cluster_true_cells)
    #print(" Cells not supported by the cluster ",cells_not_supported)
    return cluster_true_cells

''' Calculate the probability that cells exist in a cluster. '''
def calculate_prob(cells, cluster, cell_cluster_prob, cell_genotype, cluster_genotype, alpha, beta):
    cluster_genotype_list = cluster_genotype.split(',')
    #print(cluster_genotype_list)
    alpha = float(alpha)
    beta = float(beta)
    #print(" No. of cells ",len(cells))
    for cell in cells:
        p_Di_Ck = 1
        cg = cell_genotype[int(cell)].split(',')
        #print(" CG ",cg)
        key = str(cell)+'_'+str(cluster)
        for i in range(len(cluster_genotype_list)):
            #print(" Index ",i)
            #print(" Cell genotype len ",len(cg[i]))
            #print(" Cluster genotype len ",len(cluster_genotype_list[i]))
            Dij = int(cg[i])
            if Dij == 3:
                continue
            Ckj = int(cluster_genotype_list[i])
            p_Dij_Ckj = ((((1-beta)**Dij) * (beta**(1-Dij)))**Ckj) * (((1-alpha)**(1-Dij)) * (alpha ** Dij))**(1-Ckj)
            p_Di_Ck = p_Di_Ck * p_Dij_Ckj
            #print(" Probability of cell ",p_Di_Ck)
        #print(" Probability of cell not in this cluster ", (1-p_Di_Ck))
        #L_ik = np.log(p_Di_Ck / (1-p_Di_Ck))
        #print(" ======================= ")
        cell_cluster_prob[key] = p_Di_Ck
        #print(L_ik)
    #print(" ======================= ")
    return cell_cluster_prob

''' Get the wrong clones at each timepoint. '''
def get_WrongClones(cell_genotype, cluster_genotype, tp_cells, cluster_cell_tp, alpha, beta):
    tp_clusters = {}
    tp_cell_cluster_prob = {}
    tp_wrongClones = {}
    tp_true_clusters = SortedDict()
    for cluster, gen in cluster_genotype.items():
        clusterID = (cluster.split('_'))[0]
        timepoint = (cluster.split('_'))[1]
        if timepoint in tp_clusters:
            t_list = tp_clusters[timepoint]
            t_list.append(clusterID)
            tp_clusters[timepoint] = t_list
        else:
            t_list = []
            t_list.append(clusterID)
            tp_clusters[timepoint] = t_list
    #print("Timepoint clusters ",tp_clusters)
    for tp, clusters in tp_clusters.items(): # The clusters are stratified at each timepoint here.
        cell_cluster_prob = {}
        cells = tp_cells[tp]
        #print(" Cells ",cells)
        if len(clusters) == 1: # If at any timepoint there is only one clone then avoid checking likelihood.
            #print(" TP 3 cluster ",clusters[0])
            cc_dict = {}
            cc_dict[clusters[0]] = cells
            tp_true_clusters[int(tp)] = cc_dict
        else:
            for cluster in clusters:
                cell_key = cluster+'_'+tp
                cluster_cells = cluster_cell_tp[cell_key]
                #print(" Cluster ",cluster," cells",cluster_cells)
                cell_cluster_prob = calculate_prob(cells, cluster, cell_cluster_prob, cell_genotype, cluster_genotype[cell_key], alpha, beta) # Prob of each cells in the cluster.
            # print(" Timepoint ",tp," Cluster ",cell_cluster_prob.keys())
            cluster_true_cells = calculate_likelihood(cell_cluster_prob)
            tp_true_clusters[int(tp)] = cluster_true_cells
        tp_cell_cluster_prob[tp] = cell_cluster_prob
        #for cluster_true, true_cells in cluster_true_cells.items():
        #    print(cluster_true," ",len(true_cells))
        #W_clones = subclone_truly_exists(cluster_true_cells, len(cells))
        #tp_wrongClones[tp] = W_clones # Wrong clones at each timepoint.  
    #print(" Timepoint true clusters ",tp_true_clusters)
    return tp_clusters, tp_wrongClones, tp_true_clusters, tp_cell_cluster_prob

''' Reassign all remaining cells that were missed by assigning with likelihood. '''
def assign_remainingCells(tp_true_clusters, tp_cell_cluster_prob):
    #print(" Cell cluster probability ")
    #print(tp_cell_cluster_prob.keys())
    for tp, cell_cluster in tp_cell_cluster_prob.items():
        cells_list = set()
        for cells, prob in cell_cluster.items(): # Get the cells at a timepoint.
            cell = (cells.split("_"))[0]
            cells_list.add(cell)
        #print(" Timepoint ",tp," cells ",len(cells_list))

        tp_clusters = tp_true_clusters[int(tp)] 
        tp_cells_list = []
        cluster_list = set(tp_clusters.keys())
        for cluster in tp_clusters:
            tp_cells_list.extend(tp_clusters[cluster]) # Get the assigned cells.
        tp_cells_list = [int(i) for i in tp_cells_list]
        cells_list = [int(i) for i in cells_list]
        missing_cells = set(cells_list)-set(tp_cells_list) # Check the cells missed.
        if missing_cells == {}: # No need to check further if there are no missing cells at this timepoint.
            return tp_true_clusters

        print(" Cells absent at timepoint ",tp," ",missing_cells)
        for mc in missing_cells:
            missing_cell_prob = {}
            for cluster in cluster_list:
                missing_cell_prob[cluster] = cell_cluster[str(mc)+"_"+cluster]
                
            print(" Missing cell prob ",missing_cell_prob)
            max_prob_cluster = max(missing_cell_prob, key= missing_cell_prob.get) # Choose the cluster with maximum probability. If there is a tie then select a random cluster.
            old_cluster_cell_list = tp_clusters[max_prob_cluster] # Get the previously assigned cells to this cluster.
            old_cluster_cell_list.append(str(mc)) # Add the new cells to the cluster.
            tp_clusters[max_prob_cluster] = old_cluster_cell_list # Finally update the dictionary with the cells and clusters.
            #print(" Max prob cluster ",max(missing_cell_prob, key= missing_cell_prob.get))
        tp_true_clusters[int(tp)] = tp_clusters

    #print(" After reassigning cells ")
    #for tp, cluster_cell in tp_true_clusters.items():
    #    print(" Timepoint ",tp)
    #    for cluster, cell in cluster_cell.items():
    #        print(" Cluster ",cluster," # of cells ",len(cell))
    return tp_true_clusters

#''' Merge clones with correct ones based on distance. '''
#def mergeClones(tp_clusters, tp_wrongClones, cluster_cells_tp, cells_genotype, cluster_genotype, alpha, beta):
#    for tp, wrongClones in tp_wrongClones.items():
#        tp_cluster = tp_clusters[tp]
#        for wc in wrongClones:
#            cell_cluster_prob = {}
#            for cluster in tp_cluster:
#                if cluster in wrongClones: # checking only with correct clones.
#                    continue
                #correct_clone_genotype = cluster_genotype[cluster+'_'+tp].split(',')
#                wc_cells = cluster_cells_tp[wc+'_'+tp]
#                #print(" Wrong cells ",wc_cells)
#                cell_cluster_prob = calculate_prob(wc_cells, cluster, cell_cluster_prob, cells_genotype, cluster_genotype[cluster+'_'+tp], alpha, beta) # Prob of each cells in the cluster.
#            print("cell_cluster_prob  ",cell_cluster_prob)
#            closest_cluster = min(cell_cluster_prob, key=cell_cluster_prob.get) # The cluster with the smallest distance to the cells in wrong clone.
#            print(" Closest cluster ",closest_cluster)

''' After correcting the clones with probability and likelihood merge the clusters having same mutations at a timepoint. '''
def merge_tp_clusters(correced_tp_clusters, tp_true_clusters):
    #copy_tp_true_clusters = tp_true_clusters.copy()
    #print(tp_true_clusters[2]['4'])
    for tp, cluster_genotype in correced_tp_clusters.items():
        print(" Timepoint ",tp)
        cluster_to_drop = []
        for cluster1, genotype1 in cluster_genotype.items():
            g1mut_arr = [] # Get the mutations from the genotype.
            for i in range(len(genotype1)):
                if int(genotype1[i]) == 1:
                    g1mut_arr.append(i)
            for cluster2, genotype2 in cluster_genotype.items():
                if cluster1 == cluster2 or cluster1 in cluster_to_drop or cluster2 in cluster_to_drop:
                    continue
                g2mut_arr = []
                for i in range(len(genotype2)):
                    if int(genotype2[i]) == 1:
                        g2mut_arr.append(i)

                if set(g1mut_arr) == set(g2mut_arr): # Check if both clusters have same mutations
                    cluster_to_drop.append(cluster2)
                    #print(" Cluster1 ",cluster1," Cluster2 ",cluster2," Genotype ",genotype1)
                    #print(" Genotype2 ",genotype2)
                    #print("=====================================")
                    print(tp_true_clusters)
                    print(" Cluster1 ",cluster1)
                    temp_list = []
                    cluster1_cells = tp_true_clusters[tp][cluster1]
                    cluster2_cells = tp_true_clusters[tp][cluster2]
                    #print(" Cluster 1 cells ",cluster1_cells)
                    #print(" Cluster 2 cells ",cluster2_cells)
                    temp_list.extend(cluster1_cells)
                    temp_list.extend(cluster2_cells)
                    tp_true_clusters[tp][cluster1] = temp_list
                    #del copy_tp_true_clusters[tp][cluster2]
        
        for cc in cluster_to_drop: # After merging the clones and their cells drop their corresponding references.
            print(" Cluster to drop ",cc)
            del cluster_genotype[cc]
            del tp_true_clusters[tp][cc]

        #print(" After merging timepoint clusters ---- ")
        #print(cluster_genotype)
        #print(tp_true_clusters)
    return correced_tp_clusters, tp_true_clusters

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
    cell_tp = {} # Map cells to timepoints
    tp_cell = {} # Map timepoints to cells
    while(line != ""):
        line_a = line.split('\t')
        #print(line_a)
        tp = (line_a[0].split('_'))[0]
        cells = line_a[1].split(';')
        for cell in cells:
            cell_tp[cell] = tp
            if tp in tp_cell:
                c_list = tp_cell[tp]
                c_list.append(cell)
                tp_cell[tp] = c_list
            else:
                c_list = []
                c_list.append(cell)
                tp_cell[tp] = c_list
        line = file.readline().rstrip('\n')
    return cell_tp, tp_cell

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
    print(" No. of Cells in each cluster ",cc_len_time)
    print(" Cells in each cluster ",cc_time)
    return cc_time

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input Gprime matrix (Input to the Longitudinal tree)")
parser.add_argument("-cg", "--cg",dest ="cg", help="Input GDoublePrime matrix")
parser.add_argument("-D", "--D",dest ="D", help="D matrix representing true cell genotypes.")
parser.add_argument("-cc", "--cc",dest ="cc", help="Assignment file indicating cells belonging to clusters.")
parser.add_argument("-celltp", "--celltp",dest = "celltp", help="Cell timepoints.")
#parser.add_argument("-alpha","--alpha",dest = "alpha", help="False positive rate")
#parser.add_argument("-beta","--beta",dest = "beta", help="False negative rate")
parser.add_argument("-FPFN","--FPFN",dest = "FPFN", help="Estimated False positive and False negative rate from BnpC.")
parser.add_argument("-sim", "--sim",dest = "sim", help="True if comparing with simulated data else false.")
parser.add_argument("-op_tree","--op_tree",dest="op_tree", help="The output file to save the final longitudinal tree.")
parser.add_argument("-cloneNode","--cloneNode",dest="cloneNode", help="Save the numpy file mapping the clusters to the nodes in the tree.")
parser.add_argument("-cloneCells","--cloneCells",dest="cloneCells", help="Cells belonging to each cluster.")
args = parser.parse_args()
#Graph, node_tp, clone_node = inputToGraph(args.input) # Creates a graph from the input G' matrix.
#print(" Clone to nodes ",clone_node)

if args.sim == "true":
    cluster_cell = sim_cell_cluster(args.cc)
else:
    #cluster_cell = get_cell_cluster(args.cc)
    cluster_cell = np.load(args.cc,allow_pickle='TRUE').item()
#print(" Cluster cell ")
#print(cluster_cell)
#print(" ============ ")
#print(" Cluster mutation ")
cluster_mutation = get_cluster_mutations(args.input)
#print(cluster_mutation)
#print(" ============ ")
cell_tp, tp_cell = get_cell_timepoints(args.celltp)
cluster_cells_tp = getCellsOfEachCluster(cell_tp, cluster_cell)
#print(" Cell timepoints ",tp_cell)

# Correcting clusters starts here. 
# ============================================
cells_genotype = getGenotype_dict(args.D, "D")
cells_G_genotype =  getGenotype_dict(args.cg, "G") # Get the cell genotype from the consensus genotype inferred by BnpC.
#print(" Cells genotype ")
#print(cells_G_genotype)

cluster_genotype = getGenotype_dict(args.input, "cluster")
#print(" Cluster genotype ")
#print(cluster_genotype)

FPFN_dict = np.load(args.FPFN,allow_pickle='TRUE').item() # Read the estimated FP FN rate from BnpC
FP_rate = FPFN_dict['FP']
FN_rate = FPFN_dict['FN']
print("FP rate ",FP_rate,"FN rate ",FN_rate)

tp_clusters, tp_wrongClones, tp_true_clusters, tp_cell_cluster_prob = get_WrongClones(cells_genotype, cluster_genotype, tp_cell, cluster_cells_tp, FP_rate, FN_rate)
#print(" Clusters to cell dict ",cluster_cell)
#print(" Wrong clones ",tp_wrongClones)
#print(" True clusters ",tp_true_clusters)
#print(" Cell cluster prob ",tp_cell_cluster_prob)
print(" Before reassigning cells ")
for tp, cluster_cell in tp_true_clusters.items():
    print(" Timepoint ",tp)
    for cluster, cell in cluster_cell.items():
        print(" Cluster ",cluster," # of cells ",len(cell))

tp_true_clusters = assign_remainingCells(tp_true_clusters, tp_cell_cluster_prob)
correced_tp_clusters = getCorrectedClusterGenotype(cells_G_genotype, tp_true_clusters)
#print(" Corrected cluster genotype ",correced_tp_clusters)
# Merge those clusters having same mutations or consensus genotype.
correced_tp_clusters, tp_true_clusters = merge_tp_clusters(correced_tp_clusters, tp_true_clusters)

tp_cluster_cells = {}
print(" After reassigning and merging cells ")
for tp, cluster_cell in tp_true_clusters.items():
    print(" Timepoint ",tp)
    for cluster, cell in cluster_cell.items():
        print(" Cluster ",cluster," # of cells ",len(cell))
        tp_cluster_cells[cluster+"_"+str(tp)] = len(cell)

print(" Tp cluster cell dict ",tp_cluster_cells)
print(" New Graph is here ")
Graph, node_tp, clone_node = build_graph(correced_tp_clusters)
print(" Clone to nodes ",clone_node)
#mergeClones(tp_clusters, tp_wrongClones, cluster_cells_tp, cells_genotype, cluster_genotype, args.alpha, args.beta)
# =============================================
# Correcting clusters ends here.

# Longitudinal algorithm to infer the tree after clustering starts here.
up_Graph, usc_children = find_unobservedSubclones(Graph, node_tp) # Updates the graph by including unobserved subclones, adding edges and branch length.
node_s,s_pIDs = find_nodeS(usc_children)
up_Graph, usc_nodes = enroll_usc(Graph, node_s, s_pIDs)
print(" Node TP ",node_tp)
print(" Unobserved subclones at timepoints ",usc_nodes)
# Traverse the updated graph and find unconnected nodes.
up_Graph = find_nodeV(Graph, node_tp)
up_Graph = connect_danglingNodes(Graph, node_tp, usc_nodes)

for i in range(len(up_Graph)):
    print(up_Graph[i].id," ",up_Graph[i].timepoint_," ",up_Graph[i].parentID," ",up_Graph[i].type," ",up_Graph[i].mutations," ",up_Graph[i].edge_length)

#save_graph(up_Graph,"default_mr_op/rep1_tree.csv")
save_graph(up_Graph,args.op_tree)
np.save(args.cloneNode, clone_node) # Save the clone to nodes mapping to be used for evaluation
np.save(args.cloneCells, tp_cluster_cells) # Save the no of cells in each clones
