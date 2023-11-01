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
    ignore_mutations = []
    #total0s = 0
    #total1s = 0
    for j in range(0,cg_len):
        total0s = 0
        total1s = 0
        #print(" Vote cluster mutation ",j)
        for i in range(0,cell_gs_len):
            if cell_gs_list[i][j] == '3':
                #print(" For cell ",i," Ignored mutation ",j)
                ignore_mutations.append(j)
                continue
            if cell_gs_list[i][j] == '0':
                total0s = total0s+1
            if cell_gs_list[i][j] == '1':
                total1s = total1s+1
        #print(" Total 0s ",total0s," Total 1s ",total1s)
        if total0s > total1s:
            voted_cg.append('0')
        if total1s >= total0s:
            voted_cg.append('1')
    #print(" Voted genotype ",voted_cg)
    #print(" Ignore mutations ",ignore_mutations)
    return voted_cg, ignore_mutations

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

''' Check if any node in time t has same set of mutations as any of the unobserved subclones 
and eliminate those unobserved subclones. Also, eliminate unobserved subclones having 1 child node. '''
def eliminate_usc_sets(Graph, timepoint, potential_usc_sets, hash_table):
    print("====== Inside eliminating unobserved subclones ==========")
    eliminate_usc = []
    updated_potential_usc_sets = []
    # If any merged mutation has only 1 node as a value then eliminate them
    for nodes, mut in hash_table.items():
        nodes_list = nodes.split(',')
        if len(nodes_list) == 1:
            eliminate_usc.append(mut)

    # Check if any parent in time t has same set of mutations as unobserved subclones.
    # If yes then remove them.
    for i in range(len(Graph)):
        if not float(Graph[i].timepoint_) == float(timepoint):
            continue
        print(" Node id ",i)
        t_mut = Graph[i].mutations
        print(" Time t mutations ",t_mut)
        for mut in potential_usc_sets:
            print(" Mutations in USC ",mut)
            #if set(mut) < set(t_mut): # checking superset of mutations in time t
            if sorted(mut) == sorted(t_mut): # if mutation same as any node in time t
                if mut not in eliminate_usc:
                    eliminate_usc.append(mut)

    print(" Eliminated unobserved subclones ",eliminate_usc)

    for usc in potential_usc_sets:
        if usc not in eliminate_usc:
            updated_potential_usc_sets.append(usc)
    return updated_potential_usc_sets

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

''' Choose usc, parent and update the candidate unobserved subclones. '''
def filter_unobservedSubclones(Graph, time_usc, time_usc_keys, time_usc_values, t1_nodes, unlinked_usc):
    print(" ============================================= ")
    print(" Unobserved tp ",time_usc_keys)
    #max_child = max([len(time_usc[mut]) for mut in time_usc])
    #print(" Maximum children ",max_child)
    #for mut in time_usc:
    #    if len(time_usc[mut]) == max_child:
    #        usc_withMax_children = mut
    
    usc_children_list = time_usc_values
    selected_usc_mut = time_usc_keys.split(',') # The selected usc
    parent_node = "No parent"
    unlinked_usc = {}
    #parent_nodes = {} # From here we will choose the one with maximum no. of mutation intersections
    # Before selecting any parent node check if node in time k has same set of mutations as the merged mutations.
    # If yes then remove that merged mutation and don't enroll any parent node.  
    # Else select the parent node and proceed.
    for n1 in t1_nodes:
        n1_mut = Graph[n1].mutations
        #if sorted(n1_mut) == sorted(selected_usc_mut):
        #    del(time_usc[usc_withMax_children])
        #    print(n1," same as usc, hence eliminated ")
        #    return "No parent", [], [], time_usc

        print(" N1 mut ",n1_mut)
        selected_usc_mut = [int(i) for i in selected_usc_mut]
        print(" USC mut ",selected_usc_mut)
        if set(selected_usc_mut).issuperset(set(n1_mut)):
            parent_node = n1

        #intersections = set(n1_mut).intersection(set(selected_usc_mut))
        #parent_nodes[n1] = len(intersections)
        #print(n1," Intersections ",intersections)

    #max_intersection = max(parent_nodes)
    #if max_intersection == 0:
    #    parent_node = "No parent"
    #else:
    #parent_node = max(parent_nodes, key=parent_nodes.get)

    print(" Selected parent ",parent_node)
    print(" Selected children ",usc_children_list)
    
    # If there are unlinked unobserved subclones then save these before deleting them for iteration
    if parent_node == "No parent":
        unlinked_usc[time_usc_keys] = time_usc_values

    # After selecting parent, children update tp_usc to remove the selected children and unobserved subclone
    del(time_usc[time_usc_keys])
    for mut, nodes in time_usc.items():
        for n in nodes:
            if n in usc_children_list:
                nodes.remove(n)
    print(" Updated tp usc ",time_usc)
    print(" ============================================== ")
    return parent_node, usc_children_list, selected_usc_mut, time_usc, unlinked_usc

''' Check if the unobserved subclones have anymore children left. If not then delete them. '''
def recheck_unobservedSubclones(timepoint_usc):
    remove_nodes = []
    for mut, nodes in timepoint_usc.items():
        # If no children or no. of of child nodes is 1 then remove that node
        if nodes == [] or len(nodes) == 1:
            remove_nodes.append(mut)

    for usc in remove_nodes:
        del(timepoint_usc[usc])
    print(" After rechecking unobserved subclones ",timepoint_usc)
    return timepoint_usc

''' Sort the remaining unobserved subclones in decreasing order. '''
def sort_tp_usc(tp_usc):
    noOfnodes = {}
    updated_tp_usc = {}
    for mut, nodes in tp_usc.items():
        noOfnodes[mut] = len(nodes)

    sorted_noOfnodes = sorted(noOfnodes.items(), key=lambda x:x[1], reverse=True)
    #print(" Sorted usc nodes ",sorted_noOfnodes)
    ordered_usc_keys = list((dict(sorted_noOfnodes)).keys())
    #print(" Ordered usc keys ",ordered_usc_keys)
    for keys in ordered_usc_keys:
        updated_tp_usc[keys] = tp_usc[keys]
    print(" Sorted usc nodes ",updated_tp_usc)
    return updated_tp_usc

''' Calculate overlapped mutations. '''
def getOverlappedMutations(Graph, t1_nodes, unlinked_usc):
    unlinked_usc_parent = {}
    unlinked_usc_overlappedMut = {}
    for mut, nodes in unlinked_usc.items():
        usc_mut = mut.split(',')
        parent_nodes = {}
        for n1 in t1_nodes:
            n1_mut = Graph[n1].mutations
            intersections = set(n1_mut).intersection(set(usc_mut))
            parent_nodes[n1] = len(intersections)

        max_intersection = max(parent_nodes)
        if max_intersection != 0:
            parent_node = max(parent_nodes, key=parent_nodes.get)
            unlinked_usc_parent[mut] = parent_node
            unlinked_usc_overlappedMut[mut] = max_intersection

    print(" Unlinked usc parents ",unlinked_usc_parent)
    # sort the unlinked_usc_overlappedMut
    sorted_usc_overlappedMut = sorted(unlinked_usc_overlappedMut.items(), key=lambda x:x[1], reverse=True)
    print(" Unlinked usc overlapped mutations",dict(sorted_usc_overlappedMut))
    return unlinked_usc_parent, dict(sorted_usc_overlappedMut)

''' Re-check for unobserved subclones not linked. '''
def check_unlinked_usc(up_Graph, time1, t1_nodes, unlinked_usc):
    # First check the mutation of usc with parent nodes and calculate the overlapped mutations.
    # Maintain a dictionary with mutations as keys and no. of overlapped mutations as values. Sort this dict decreasingly.
    # Loop over the sorted dict to find the parent in time k with max. no. of overlapped/intersection of mutations.
    # Remove the keys and values while looping over this in a similar way to the previous approach.
    unlinked_usc_parent, usc_overlappedMut = getOverlappedMutations(up_Graph, t1_nodes, unlinked_usc)
    for usc in usc_overlappedMut:
        if usc not in unlinked_usc:
            continue
        parent_node = unlinked_usc_parent[usc]
        usc_children = unlinked_usc[usc]
        usc_mutation = usc.split(',')
        usc_mutation = [int(i) for i in usc_mutation]
        # Remove these usc_children from other children of usc. Then if any usc has no children or just 1 child then remove that usc.
        del(unlinked_usc[usc]) # since it is already selected we remove this
        remove_usc = []
        for mut, nodes in unlinked_usc.items():
            for n in nodes:
                if n in usc_children:
                    nodes.remove(n)

            if nodes == [] or len(nodes) == 1:
                remove_usc.append(mut)

        for ru in remove_usc:
            del(unlinked_usc[ru])

        # Now add the parent and child nodes in the Graph
        usc_nodeId = len(up_Graph)
        up_Graph.append(graphnode(usc_nodeId))
        up_Graph[usc_nodeId].id = usc_nodeId

        p_mut = up_Graph[int(parent_node)].mutations
        #usc_mutation_forTree = []
        #usc_mutation_forTree.extend(p_mut)
        #usc_mutation_forTree.extend(usc_mutation)
        #print(" USC mutation for tree ",set(usc_mutation_forTree))
        #up_Graph[usc_nodeId].mutations = list(set(usc_mutation_forTree))
        up_Graph[usc_nodeId].mutations = list(set(usc_mutation))

        up_Graph[usc_nodeId].type = 'unobserved subclone'
        usc_child_timepoint = float(up_Graph[int(usc_children[0])].timepoint_)
        up_Graph[usc_nodeId].timepoint_ = (float(time1) + usc_child_timepoint) / 2 # Update usc timepoints.

        pID = []
        pID.append(int(parent_node)) # Get the parent node Id
        up_Graph[usc_nodeId].parentID = pID
        edge_len = []
        edge_len.append(len(set(usc_mutation) - set(p_mut)))
        up_Graph[usc_nodeId].edge_length = edge_len # New mutations in the usc when compared with node in t.

        for usc_child in usc_children:
            usc_child_pIds = [usc_nodeId] # update the parent node to include usc as its parent
            up_Graph[int(usc_child)].parentID = usc_child_pIds

            usc_child_mut = up_Graph[int(usc_child)].mutations
            usc_child_edge = up_Graph[int(usc_child)].edge_length
            usc_child_edge.append(len(set(usc_child_mut) - set(usc_mutation)))
            up_Graph[int(usc_child)].edge_length = usc_child_edge # No. of new mutations in node t+1 compared to usc.

    return up_Graph

''' Select the unobserved subclone which has the largest number of child nodes.
To connect with the parent we choose the node in time t that has largest intersection of mutations. '''
def select_unobservedSubclones(up_Graph, node_tp, usc_children, secondLevel):
    # For second level of unobserved subclone check with a flag secondLevel.
    # If secondLevel then t1_nodes should have nodes from previous two timepoints.
    print(" Maximum timepoint ",max(node_tp))
    print(" USC children ",usc_children)
    for time1, nodes in node_tp.items():
        # we don't check the last tp because all unobserved subclones have already occured
        if time1 == max(node_tp) or time1 not in usc_children:
            continue
        print("USC Timepoint ",time1," Nodes ",nodes)
        if secondLevel:
            t1_nodes = []
            if time1 == 1:
                continue
            prev_timepoint_nodes = node_tp[time1-0.5]
            print(time1," Second level parent timepoints ",time1-0.5," nodes ",prev_timepoint_nodes)
            # Include only those parent nodes are not connected to any unobserved subclones in the first layer
            for ptn in prev_timepoint_nodes:
                ptn_parent_node = up_Graph[ptn].parentID # get parent of this node
                #if ptn_parent_node == [] or ptn_parent_node == None:
                #    continue
                #ptn_parent_node_tp = up_Graph[ptn_parent_node[0]].timepoint_ # get the timepoint of the parent
                ptn_parent_node_tp = time1-0.5
                if ptn_parent_node_tp != time1:
                    t1_nodes.extend(prev_timepoint_nodes)
            t1_nodes.extend(nodes)
            t1_nodes = list(set(t1_nodes))
            print(" Second level parent nodes ",t1_nodes)
        else:
            t1_nodes = nodes

        tp_usc = usc_children[time1]
        # Sort the timepoint usc
        if tp_usc != {}:
            tp_usc = sort_tp_usc(tp_usc)

        copy_tp_usc = copy.deepcopy(tp_usc) # Make changes in the deepcopy
        unlinked_usc = {}
        # a while loop to check if tp_usc is empty
        # have a function to choose the usc with highest children and find the parent with largest intersection. 
        # In the same function above update tp_usc by removing child nodes from other usc and removing the observed usc.
        # Then update the graph to include parent, children
        while(copy_tp_usc != {}):
            # Also check if there are some children of the unobserved children. If not then remove all the unobserved subclones as well.
            copy_tp_usc = recheck_unobservedSubclones(copy_tp_usc)
            if copy_tp_usc == {}:
                break
            mut = list(copy_tp_usc.keys())[0]
            nodes = list(copy_tp_usc.values())[0]
            print(" Selected Mutation ",mut," Nodes ",nodes)
            parent_node, usc_children_list, usc_mutation, copy_tp_usc, unlinked_usc = filter_unobservedSubclones(up_Graph, copy_tp_usc, mut, nodes, t1_nodes, unlinked_usc)

            if parent_node == "No parent" or usc_children_list == []: # We don't enroll anymore unobserved subclones if there are no parents
                continue

            # Next update the Graph with the new unobserved subclones along with their parent and children
            usc_nodeId = len(up_Graph)
            up_Graph.append(graphnode(usc_nodeId))
            up_Graph[usc_nodeId].id = usc_nodeId

            p_mut = up_Graph[int(parent_node)].mutations
            #usc_mutation_forTree = []
            #usc_mutation_forTree.extend(p_mut)
            #usc_mutation_forTree.extend(usc_mutation)
            #print(" USC mutation for tree ",set(usc_mutation_forTree))
            #up_Graph[usc_nodeId].mutations = list(set(usc_mutation_forTree))
            up_Graph[usc_nodeId].mutations = usc_mutation

            up_Graph[usc_nodeId].type = 'unobserved subclone'
            usc_child_timepoint = float(up_Graph[int(usc_children_list[0])].timepoint_)
            up_Graph[usc_nodeId].timepoint_ = (float(time1) + usc_child_timepoint) / 2 # Update usc timepoints.

            pID = []
            pID.append(int(parent_node)) # Get the parent node Id
            up_Graph[usc_nodeId].parentID = pID
            edge_len = []
            #edge_len.append(len(set(usc_mutation_forTree) - set(p_mut)))
            edge_len.append(len(set(usc_mutation) - set(p_mut)))
            up_Graph[usc_nodeId].edge_length = edge_len # New mutations in the usc when compared with node in t.

            for usc_child in usc_children_list:
                usc_child_pIds = [usc_nodeId] # update the parent node to include usc as its parent
                up_Graph[int(usc_child)].parentID = usc_child_pIds

                usc_child_mut = up_Graph[int(usc_child)].mutations
                usc_child_edge = up_Graph[int(usc_child)].edge_length
                #usc_child_edge.append(len(set(usc_child_mut) - set(usc_mutation_forTree)))
                usc_child_edge.append(len(set(usc_child_mut) - set(usc_mutation)))
                up_Graph[int(usc_child)].edge_length = usc_child_edge # No. of new mutations in node t+1 compared to usc.

        # Checking for unlinked unobserved subclonies
        if unlinked_usc != {}:
            up_Graph = check_unlinked_usc(up_Graph, time1, t1_nodes, unlinked_usc)

    print(" Graph after updating unobserved subclones ================== ")
    print(" Also update the node_tp for second level of unobserved node check ")
    for i in range(len(up_Graph)):
        print(up_Graph[i].id," ",up_Graph[i].timepoint_," ",up_Graph[i].parentID," ",up_Graph[i].type," ",up_Graph[i].mutations," ",up_Graph[i].edge_length)
        if up_Graph[i].timepoint_ in node_tp:
            node_list = node_tp[up_Graph[i].timepoint_]
            if up_Graph[i].id not in node_list:
                node_list.append(up_Graph[i].id)
            node_tp[up_Graph[i].timepoint_] = node_list
        else:
            node_list = [up_Graph[i].id]
            node_tp[up_Graph[i].timepoint_] = node_list    

    print(node_tp)
    return up_Graph, node_tp

''' Build the hash table with nodes and mutations to reduce complexity while listing the subsets. '''
def build_hash_table(all_n2_mutations):
    mut_dict = {} # Hash table 1 where mutation is the key and nodes are values
    hash_table = {} # Hash table 2 which have merged mutations with same set of nodes. Nodes are keys and mutations values.
    #merge_mut = []
    for node, mut_list in all_n2_mutations.items():
        for mut in mut_list:
            if mut in mut_dict:
                temp_list = mut_dict[mut]
                temp_list.append(node)
                mut_dict[mut] = temp_list
            else:
                mut_dict[mut] = [node]
    
    print(" Inside build hash table mut dict ",mut_dict)
    print(" All n2 mutations ",all_n2_mutations)
    # Merge the mutations having same set of nodes
    #for mut, nodes in mut_dict.items():
    #    nodes.sort() # sorting before comparing the nodes
    #    merged_mut_list = [mut]
    #    for mut1, nodes1 in mut_dict.items():
    #        if mut == mut1:
    #            continue
    #        nodes1.sort()
    #        if nodes == nodes1:
    #            if mut1 not in merged_mut_list:
    #                merged_mut_list.append(mut1)
    #    merged_mut_list.sort()
    #    if merged_mut_list not in merge_mut:
    #        merge_mut.append(merged_mut_list)
    #print(" Merged mutations ",merge_mut)
    
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
 
''' From the final set of unobserved subclones we already know their children.Just update the hash table to exclude the eliminated unobserved subclones. '''
def get_usc_children(hash_table, final_possible_usc):
    usc_children = {}
    for nodes, mut in hash_table.items():
        if mut in final_possible_usc:
            mut_key = ",".join(map(str,mut))
            usc_children[mut_key] = nodes.split(',')
    return usc_children

''' Merge mutations for unobserved subclones to get the updated mutations. '''
def merge_mutations(hash_table):
    for nodes, mut in hash_table.items():
        node = nodes.split(',')
        for nodes1, mut1 in hash_table.items():
            print(" Nodes ",nodes1," mutations ",mut1)
            if nodes == nodes1:
                continue
            node1 = nodes1.split(',')
            if set(node1).issubset(set(node)):
                merged_mut = []
                merged_mut.extend(mut1)
                merged_mut.extend(mut)
                hash_table[nodes1] = list(set(merged_mut)) 

    print(" Updated hash table ",hash_table)
    return hash_table

''' Finding unobserved subclones by building a hash table. '''
def find_unobservedSubclones(Graph, node_tp, secondLevel):
    # For second level of checking unobserved subclones use a boolean flag secondLevel.
    # If secondLevel then if time1 an integer value then continue. time2 should be +0.5.
    # Skip the step to connect nodes with similar nodes for secondLevel.
    tp_usc_children = {} # Update unobserved subclone for each timepoint
    for time1, nodes in node_tp.items():
        if secondLevel:
            if float(time1).is_integer(): # Here we only check for unobserved subclones.
                continue
            time2 = time1+0.5
        else:
            time2 = time1+1
        t1_nodes = nodes
        if time2 in node_tp: # Compare nodes in t with nodes in t+1.
            t2_nodes = node_tp[time2]
        else:
            continue

        print(" Timepoint ==== ",time1)
        #usc_timepoint = float((time1+(time1+1))/2) # Timepoint for the unobserved subclones.
        #print(" Usc timepoint ",usc_timepoint)
        connected_n2 = []
        #all_n2_mutations = []
        all_n2_mutations = {} # use this to build the hash table 

        for n1 in t1_nodes:
            n1_mut = Graph[n1].mutations
            #print(" Node 1 ",n1," N1 mut ",n1_mut)
            #start = time.time()
            
            for n2 in t2_nodes: # Compare mutations of nodes between two timepoints to find unobserved subclones.
                if n2 in connected_n2:
                    continue
                n2_mut = Graph[n2].mutations
                if n1_mut == n2_mut:
                    if not secondLevel:
                        pID_list = []
                        pID_list.append(n1)
                        Graph[n2].parentID = pID_list
                        Graph[n2].edge_length = [0] # since there are no new mutations, the number is 0.

                        connected_n2.append(n2) # Mark connected nodes in t and t+1.
                    if n2 in all_n2_mutations: # Second check to see if same nodes exist in t and t+1 then don't find their supersets.
                        del(all_n2_mutations[n2])
                    continue

                if n2 not in all_n2_mutations:
                    all_n2_mutations[n2] = n2_mut
                    #all_n2_mutations.append(n2_mut)
            #print(" Length of Graph ",len(Graph))

        # Build the hash table here. This is for one layer of unobserved subclones.
        hash_table = build_hash_table(all_n2_mutations)
        print(" Before merging hash table ",hash_table)
        hash_table = merge_mutations(hash_table)
        print(" After merging hash table ",hash_table)
        potential_usc_sets = list(hash_table.values())
        final_possible_usc = eliminate_usc_sets(Graph, time1, potential_usc_sets, hash_table)

        # Next we eliminate those unobserved sets that are subset of any node in time k as it will not contain any new mutation.
        #eliminate_usc_sets(Graph, time1, potential_usc_sets)

        #if time1 == 1:
        #    final_possible_usc = potential_usc_sets # At t=1 or when root we directly use the sets we get without looking for supersets. 
        #    print(" Final possible usc ",final_possible_usc)
        #else:
            #final_possible_usc = find_supersets_tNodes(Graph, time1, potential_usc_sets) 
        #    final_possible_usc = eliminate_usc_sets(Graph, time1, potential_usc_sets, hash_table)
        print(" Final possible usc ",final_possible_usc)

        #usc_children = get_usc_children(Graph, time1, final_possible_usc, usc_children) # Update the usc node.
        usc_children = get_usc_children(hash_table, final_possible_usc)
        print(" USC children ",usc_children)
        tp_usc_children[time1] = usc_children

    #print(" Timepoint usc dict ",len(timepoint_usc[2]))
    return Graph, tp_usc_children

''' Look for node v such that the mutations of v that are not 
    existent in u is the least, whereas v is a node in
    time point k, or the unobserved nodes between time points k and k + 1.'''
def select_nodeV_withLeastMutation(v_mut, u_mut):
    v_mut_diff = {}
    for v, mut in v_mut.items():
        diff_mut = set(mut) - set(u_mut)
        v_mut_diff[v] = len(diff_mut)

    nodeV_withMinMutationsDiff = min(v_mut_diff,key=v_mut_diff.get)
    return nodeV_withMinMutationsDiff

''' For such nodes u we find a node v in time k in the following way:
    1. For such nodes u we find a node v in time k whose set of mutations is a subset
of those in u. 
    2. If such a node v does not exist, we connect u with node v such that 
    the mutations of v that are not existent in u is the least, whereas v is a node in
time point k, or the unobserved nodes between time points k and k + 1.
If there are more than one node v then we select one randomly. '''
def find_nodeV(Graph, node_tp):
    timepoints = node_tp.keys()
    print(" ============== Inside finding node v ============ ")
    for i in range(1,len(timepoints)):
        t2 = timepoints[i]
        if not float(t2).is_integer(): # we only check for nodes in time k+1 or integer timepoints
            continue
        t1 = t2 - 1
        usc_nodes_tp = t2 - 0.5

        print(" Time ",timepoints[i])
        t2_nodes = node_tp[t2]
        t1_nodes = node_tp[t1]
        if usc_nodes_tp in node_tp:
            usc_nodes = node_tp[usc_nodes_tp]
        else:
            usc_nodes = []

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
            
            # If v node is not existent then we do the following steps.
            if v_nodes == []:
                t1_nodes.extend(usc_nodes) # We also check for unobserved subclones hence add them here.
                v_mut_dict = {}
                for v in t1_nodes:
                    v_mut = Graph[v].mutations
                    v_mut_dict[v] = v_mut
                    
                selected_nodeV = select_nodeV_withLeastMutation(v_mut_dict, u_mut)
                v_nodes.append(selected_nodeV)
                v_mutations.append(v_mut_dict[selected_nodeV])
                #largest_v_mut = find_mostConnectedNodes(v_mutations) # Find the largest mutation set size.
                #print(" largest_v_mut ",largest_v_mut)
            #if largest_v_mut == []:
            #    continue

            #if len(largest_v_mut) > 1:
            #    random_index = random.randrange(len(largest_v_mut)) # Find a random mutation set if more than one exist.
            #    largest_v_mut = largest_v_mut[random_index]
            #else:
            #    largest_v_mut = largest_v_mut[0]

            #for i in range(len(v_mutations)):
            #    if set(largest_v_mut) == set(v_mutations[i]):
            #        v_node_index = i
            #        break

            pID = []
            pID.extend(v_nodes)
            Graph[u].parentID = pID

            edge_len = []
            edge_len.append(len(set(u_mut) - set(v_mutations[0])))
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

''' Even after correcting and merging clones there might be few clones with just 1 cell. We merge them randomly with any one of the clones. '''
def merge_clones_WithOneCell(correced_tp_clusters, tp_true_clusters):
    cells_to_merge = {} # Maintain a dict to merge only 1 cells
    print(" Before merging 1 cells into clusters ")
    print(tp_true_clusters)

    for tp, cluster_cell in tp_true_clusters.items():
        print(" Timepoint ",tp)
        cluster_to_drop = []
        for cluster, cell in cluster_cell.items():
            print(" Cluster ",cluster," # of cells ",len(cell))
            if len(cell) == 1: # First look for the cluster with just 1 cell.
                cluster_to_drop.append(cluster)
                if tp in cells_to_merge:
                    cells_list = cells_to_merge[tp]
                    cells_list.append(cell)
                    cells_to_merge[tp] = cells_list
                else:
                    cells_to_merge[tp] = [cell]

        # Then drop that cluster
        for cc in cluster_to_drop:
            print(" Cluster to drop ",cc)
            del correced_tp_clusters[tp][cc]
            del tp_true_clusters[tp][cc]

    # Reassign the cells to the first existing cluster
    for tp in tp_true_clusters:
        if tp not in cells_to_merge:
            continue
        cells_list = cells_to_merge[tp]
        clusters = list(tp_true_clusters[tp].keys())
        print(" Clusters ",clusters)
        existingCells = tp_true_clusters[tp][clusters[0]]
        existingCells.append(cells_list)
        tp_true_clusters[tp][clusters[0]] = existingCells

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

''' Get # of cells in each cluster. Stratification of clusters in timepoint. '''
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
        #print(" Likelihood ",L_ik)
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

    print(" Clusters true cells ",cluster_true_cells)
    print(" Cells not supported by the cluster ",cells_not_supported)
    return cluster_true_cells

''' Re-determine number of subclones and consensus of subclones. 
tp_prev_unassigned_cells are the unassigned cells before re-calculating FP FN. '''
def updateSubclones(tp_true_clusters, tp_unassigned_cells, tp_prev_unassigned_cells, tp_cell_cluster_prob, D_cell_genotype):
    tp_subclone_genotype = {}
    print(" Inside updated subclones code ========= ")
    for tp, true_clusters in tp_true_clusters.items():
        unassigned_cells = tp_unassigned_cells[tp]
        prev_unassigned_cells = tp_prev_unassigned_cells[str(tp)]
        prev_unassigned_cells = [str(i) for i in prev_unassigned_cells]
        # update the unassigned_cells to include the prev_unassigned_cells
        #unassigned_cells.extend(prev_unassigned_cells)
        #if unassigned_cells == [] or prev_unassigned_cells == []:
        #    continue
        print(tp," True cluster ",true_clusters)

        cell_cluster_prob = tp_cell_cluster_prob[tp]
        subclone_genotype = {}
        skip_cluster = []
        
        # First get the subclone consensus genotype
        for cluster, cells in true_clusters.items():
            cell_genotype = []
            # Since no. of supporting cells less than 2 we eliminate this subclone
            print(" Cells ",cells)
            print(" Prev unassigned cells ",prev_unassigned_cells)
            if len(cells) < 2: 
                unassigned_cells.extend(cells)
                skip_cluster.append(cluster) # eliminate this cluster
                #del(tp_true_clusters[tp][cluster])
                continue
            for cell in cells:
                if cell in prev_unassigned_cells: # don't consider these cells while getting the voted genotype
                    print(" Cell not considered for majority vote ",cell)
                    continue
                cell_genotype_list = D_cell_genotype[int(cell)].split(',')
                cell_genotype.append(cell_genotype_list)
            #if tp == 3:
            #    print(" ============================= ")
            #    print(" Timepoint ",tp," Cluster ",cluster)
            #    print(cell_genotype)
            print(" Timepoint ",tp," Cluster ",cluster) 
            subclone_genotype[cluster], imut = vote_cluster_genotypes(cell_genotype)
        tp_subclone_genotype[tp] = subclone_genotype # Voted subclone with the supporting cells

        # Drop the eliminated clusters
        for sc in skip_cluster:
            print(tp," Dropping cluster ",sc)
            del(tp_true_clusters[tp][sc])

        # Now re-assign the unassigned cells with tie or existing alone in a subclone
        cells_list = []
        cluster_list = []
        for cc, prob in cell_cluster_prob.items():
            cell = cc.split('_')[0]
            cells_list.append(cell)
            cluster = cc.split('_')[1]
            cluster_list.append(cluster)

        print(" Unassigned cells ",unassigned_cells)
        for uc in unassigned_cells:
            for cell in set(cells_list):
                if uc != cell:
                    continue
                cell_prob = {}
                for cluster in set(cluster_list):
                    print(" Cluster ",cluster," Skip cluster ",skip_cluster," Cluster keys ",list(true_clusters.keys()))
                    print(" Unassigned cell ",uc)
                    # Don't check eliminated clusters
                    if cluster in skip_cluster or cluster not in list(true_clusters.keys()):
                        print(" Inside the If loop ")
                        continue
                    # cells probability for each cluster
                    cell_prob[cluster] = cell_cluster_prob[str(cell)+"_"+cluster]
                print(" Cell cluster prob ",cell_prob)
                #if cell_prob == {}:
                #    return "Don't check further", "Don't check further"
                # Choose the cluster with highest likelihood
                max_prob_cluster = max(cell_prob, key=cell_prob.get)
                print(" Max prob cluster ",max_prob_cluster)
                print(" True cluster keys ",true_clusters.keys())
                existing_cells = true_clusters[max_prob_cluster]
                existing_cells.append(cell) # Add unassigned cell to existing cells
                true_clusters[max_prob_cluster] = existing_cells

    print(" Updated subclone genotype ",tp_subclone_genotype)
    return tp_subclone_genotype, tp_true_clusters

''' Re-assign cells according to maximum probability or likelihood. '''
def assignCells_MaxProb(cell_cluster_prob):
    cells_list = []
    unassigned_cells_list = []
    cluster_list = []
    cluster_true_cells = {}

    for cc, prob in cell_cluster_prob.items():
        cell = cc.split('_')[0]
        cells_list.append(cell)
        cluster = cc.split('_')[1]
        cluster_list.append(cluster)

    for cell in set(cells_list):
        cell_prob = {}
        for cluster in set(cluster_list):
            # cells probability for each cluster
            cell_prob[cluster] = cell_cluster_prob[str(cell)+"_"+cluster]

        #print(cell," cell prob ",cell_prob)
        max_prob = max(cell_prob.values())
        #print(" Maximum prob ",max_prob)
        max_prob_clusters = []
        for cluster, prob in cell_prob.items():
            if prob == max_prob:
                max_prob_clusters.append(cluster)

        if len(max_prob_clusters) > 1: # This means that there is a tie between clusters
            # If there is tie then we don't assign these cells to any cluster now.
            #print(" Tied clusters ",max_prob_clusters)
            unassigned_cells_list.append(cell)
        else:
            max_prob_cluster = max(cell_prob, key=cell_prob.get)
            #print(" Max prob cluster ",max_prob_cluster)
            if max_prob_cluster in cluster_true_cells:
                temp_cells = cluster_true_cells[max_prob_cluster]
                temp_cells.append(cell)
                cluster_true_cells[max_prob_cluster] = temp_cells
            else:
                cluster_true_cells[max_prob_cluster] = [cell]
        #print(" Cell ",cell," cluster ",max_prob_cluster)
    print(" Clusters true cells ",cluster_true_cells)

    #total_cells = 0 # check that total no. of cells is intact
    #for cluster, cells in cluster_true_cells.items():
    #    total_cells = total_cells+len(cells)
    print(" Unassigned cells ",unassigned_cells_list)
    return cluster_true_cells, unassigned_cells_list

''' Calculate the probability that cells exist in a cluster. '''
def calculate_prob(cells, cluster, cell_cluster_prob, cell_genotype, cluster_genotype, ignore_mutations, alpha, beta):
    print(" Cluster genotype ",cluster_genotype)
    #cluster_genotype_list = cluster_genotype.split(',')
    #print(" Cluster mutations ",cluster_genotype_list)
    #print(" Ignore mutations ",ignore_mutations.keys())
    alpha = float(alpha)
    beta = float(beta)
    # Total no. of mutations after removing the missing mutations
    #total_mutations = len(cluster_genotype) - len(set(ignore_mutations))
    #print(" Total mutations ",total_mutations)
    #print(" No. of cells ",len(cells))
    for cell in cells:
        p_Di_Ck = 1 # i = cell, j = mutation, k = subclone
        #p_Di_Ck = 0
        cg = cell_genotype[int(cell)].split(',')

        cell_ignore_mutations = ignore_mutations[int(cell)]
        total_mutations = len(cluster_genotype) - len(set(cell_ignore_mutations))
        print(" Total mutations ",total_mutations)
        #if cell == '51' or cell == '52':
        #    print(" CG ",cg)
        #    print(" consensus genotype ",cluster_genotype_list)
        key = str(cell)+'_'+str(cluster)
        for j in range(len(cluster_genotype)):
            #print(" Index ",i)
            #print(" Cell genotype len ",len(cg[i]))
            #print(" Cluster genotype len ",len(cluster_genotype_list[i]))
            Dij = int(cg[j])
            if Dij == 3 or j in cell_ignore_mutations:
                #if i not in ignore_mutations:
                #    print(" Not matched mutations ",i)
                continue
            Ckj = int(cluster_genotype[j])
            #if cell == '51' or cell == '52':
            #    print(" For cell ============== ",cell)
            #    print(Dij)
            #    print(Ckj)
            p_Dij_Ckj = ((((1-beta)**Dij) * (beta**(1-Dij)))**Ckj) * (((1-alpha)**(1-Dij)) * (alpha ** Dij))**(1-Ckj)
            p_Di_Ck = p_Di_Ck * p_Dij_Ckj
            #p_Di_Ck = p_Di_Ck + np.log(p_Dij_Ckj)
            #print(" Probability of cell ",p_Di_Ck)
        #print(" Probability of cell not in this cluster ", (1-p_Di_Ck))
        #L_ik = np.log(p_Di_Ck / (1-p_Di_Ck))
        #print(" ======================= ")
        #cell_cluster_prob[key] = np.log(p_Di_Ck) / total_mutations
        if total_mutations != 0:
            cell_cluster_prob[key] = p_Di_Ck / total_mutations
        else:
            cell_cluster_prob[key] = 0
        #print(key," probability ",cell_cluster_prob[key])
    #print(" ======================= ")
    return cell_cluster_prob

''' Get the wrong clones at each timepoint. '''
def get_WrongClones(cell_genotype, D_cluster_genotype, tp_cells, cluster_cell_tp, ignore_mutations, tp_alpha, tp_beta, tp_prev_unassigned_cells):
    tp_clusters = {}
    tp_unassigned_cells = {}
    tp_cell_cluster_prob = {}
    tp_wrongClones = {}
    tp_true_clusters = SortedDict()

    # Our voted Dmatrix is cluster genotype here
    for cluster, gen in D_cluster_genotype.items(): 
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

    print("Timepoint clusters ",tp_clusters)
    for tp, clusters in tp_clusters.items(): # The clusters are stratified at each timepoint here.
        print(" Timepoint ",tp)
        print(" ============================================= ")
        alpha = tp_alpha[tp]
        beta = tp_beta[tp]
        cell_cluster_prob = {}
        cells = tp_cells[tp]
        prev_unassigned_cells = tp_prev_unassigned_cells[tp]
        # If there are no unassigned cells then we don't change anything
        #if prev_unassigned_cells == []:
        #    tp_true_clusters[int(tp)] = cells
        #    tp_cell_cluster_prob, 
        #    tp_unassigned_cells[int(tp)] = []
        #    continue

        #print(" Cells ",cells)
        if len(clusters) == 1: # If at any timepoint there is only one clone then avoid checking likelihood and place all cells there.
            #print(" TP 3 cluster ",clusters[0])
            cc_dict = {}
            #print(" All cells ",cells)
            #print(" Prev unassigned cells ",prev_unassigned_cells)
            #prev_unassigned_cells = [str(i) for i in prev_unassigned_cells]
            #new_cells = set(cells) - set(prev_unassigned_cells)
            #print(" New cells ",len(new_cells))
            cc_dict[clusters[0]] = cells
            tp_true_clusters[int(tp)] = cc_dict
            tp_unassigned_cells[int(tp)] = []
        else:
            for cluster in clusters:
                cell_key = cluster+'_'+tp
                cluster_cells = cluster_cell_tp[cell_key]
                print(" Inside get Wrong clones =========== ")
                print(" Cluster ",cluster," cells",cluster_cells)
                cell_cluster_prob = calculate_prob(cells, cluster, cell_cluster_prob, cell_genotype, D_cluster_genotype[cell_key], ignore_mutations, alpha, beta) # Prob of each cells in the cluster.
            print(" Timepoint ",tp," Cluster cell prob ",cell_cluster_prob.keys())
            cluster_true_cells, unassigned_cells_list = assignCells_MaxProb(cell_cluster_prob)
            #cluster_true_cells = calculate_likelihood(cell_cluster_prob)
            tp_true_clusters[int(tp)] = cluster_true_cells
            tp_unassigned_cells[int(tp)] = unassigned_cells_list
        tp_cell_cluster_prob[int(tp)] = cell_cluster_prob
        #for cluster_true, true_cells in cluster_true_cells.items():
        #    print(cluster_true," ",len(true_cells))
        #W_clones = subclone_truly_exists(cluster_true_cells, len(cells))
        #tp_wrongClones[tp] = W_clones # Wrong clones at each timepoint.
    #print(" Timepoint true clusters ",tp_true_clusters)
    return tp_clusters, tp_wrongClones, tp_true_clusters, tp_cell_cluster_prob, tp_unassigned_cells

''' Get the FP and FN by comparing our voted D Consensus genotype and all cells from D matrix. '''
def calc_FPFN(G_list, D_list):
    est_FP = 0
    est_FN = 0
    TN = 0
    TP = 0
   
    print(" D_list ",D_list)
    for i in range(len(D_list)):
        for j in range(len(D_list[0])):
            if D_list[i][j] == '3': # ignore the missing mutations
                continue

            if G_list[j] == '1' and D_list[i][j] == '0': # G matrix here is the inferred genotype from BnpC and we are comparing it with D.
                est_FN = est_FN+1
            if D_list[i][j] == '1' and G_list[j] == '0':
                est_FP = est_FP+1
            if D_list[i][j] == '0' and G_list[j] == '0':
                TN = TN+1
            if D_list[i][j] == '1' and G_list[j] == '1':
                TP = TP+1

    print(" Total estimated FP ",est_FP," True negatives ",TN)
    print(" Total estimated FN ",est_FN," True positives ",TP)
    #print(" FP rate deno ",est_FP+TP)
    #print(" FN rate deno ",est_FN+TN)

    est_FP_rate = est_FP/(est_FP+TP)
    est_FN_rate = est_FN/(est_FN+TN)
    rounded_est_FP_rate = round(est_FP_rate,5)
    rounded_est_FN_rate = round(est_FN_rate,5)
    print(" Estimated FP rate ",rounded_est_FP_rate)
    print(" Estimated FN rate ",rounded_est_FN_rate)
    return rounded_est_FP_rate, rounded_est_FN_rate

''' Get ignore mutations specific to cells. '''
def getCellsMissingMutations(D_matrix):
    cells_ignore_mutation = {}
    for i in range(len(D_matrix)):
        D_cellgen = D_matrix[i].split(',')
        ignore_mut = []
        for j in range(len(D_cellgen)):
            if D_cellgen[j] == '3':
                ignore_mut.append(j)
        cells_ignore_mutation[i] = ignore_mut
    print(" Cells ignore mut ",len(cells_ignore_mutation))
    return cells_ignore_mutation

''' Get the FP, FN rate for each timepoint based on majority voting of
consensus genotype by D. '''
def getTimepoint_FPFN(cluster_cells_tp, D_matrix):
    tp_FP = {} # Save the FP rate for each timepoint
    tp_FN = {} # Save the FN rate for each timepoint
    voted_D_cg = {} # Use this to get the consensus genotype with voted D
    #cluster_ignore_mutations = {} # Save the mutations to ignore while calculating the likelihood
    cells_ignore_mutations = getCellsMissingMutations(D_matrix)

    for cluster, cells in cluster_cells_tp.items():
        tp = cluster.split('_')[1]
        print(" For cluster ",cluster)
        all_cells_D = []
        for cell in cells:
            D_cellgen = D_matrix[cell].split(',')
            all_cells_D.append(D_cellgen)

        # Make our own consensus by majority vote on D
        voted_subclone_genotype, ignore_mut = vote_cluster_genotypes(all_cells_D)
        voted_D_cg[cluster] = voted_subclone_genotype
        print(" Voted D CG ",voted_subclone_genotype)
        #cluster_ignore_mutations[cluster] = ignore_mut

        if len(cells) >= 3:
            fp_rate, fn_rate = calc_FPFN(voted_subclone_genotype, all_cells_D)

            if tp in tp_FP or tp in tp_FN:
                fp_list = tp_FP[tp]
                fn_list = tp_FN[tp]
                fp_list.append(str(fp_rate)+'_'+str(len(cells)))
                fn_list.append(str(fn_rate)+'_'+str(len(cells)))
                tp_FP[tp] = fp_list
                tp_FN[tp] = fn_list
            else:
                tp_FP[tp] = [str(fp_rate)+'_'+str(len(cells))]
                tp_FN[tp] = [str(fn_rate)+'_'+str(len(cells))]

    print(" Timepoint FPs ",tp_FP)
    print(" Timepoint FNs ",tp_FN)
    # Calculate the weighted FP anf FNs for each timepoints.
    tp_alpha = {}
    tp_beta = {}
    for tp, fps in tp_FP.items():
        print(" Timepoint ",tp)
        fns_list = tp_FN[tp]
        wght_avg_FPnum = 0
        wght_avg_FNnum = 0
        wght_avg_deno = 0

        for i in range(len(fps)):
            fp_rate = float(fps[i].split('_')[0])
            fn_rate = float(fns_list[i].split('_')[0])
            cells = int(fps[i].split('_')[1])
            wght_avg_FPnum = wght_avg_FPnum + (cells * fp_rate)
            wght_avg_FNnum = wght_avg_FNnum + (cells * fn_rate)
            wght_avg_deno = wght_avg_deno + cells
            print(" FP rate ",fp_rate," FN rate ",fn_rate," cells ",cells)

        weighted_FP = wght_avg_FPnum / wght_avg_deno
        weighted_FN = wght_avg_FNnum / wght_avg_deno
        tp_alpha[tp] = weighted_FP
        tp_beta[tp] = weighted_FN
        print(" Timepoint ",tp," FP ",weighted_FP," FN ",weighted_FN)
    
    print(" Voted CG ",voted_D_cg)
    #print(" cluster ignore mut ",cluster_ignore_mutations)
    return tp_alpha, tp_beta, voted_D_cg, cells_ignore_mutations

''' Get the missing rate for each timepoint based on no. of 3s /  total number of cells * mutations in that timepoint. '''
def timepoint_missingRate(tp_cell, D_matrix):
    tp_D = {} # Get all the entries for each timepoint
    tp_MR = {} # Missing rate at each timepoint
    #print(" Timepoint missing rate ")

    for tp, cells in tp_cell.items():
        #print(" timepoint ",tp,"cells ",cells)
        all_cells_D = []
        for cell in cells:
            D_cellgen = D_matrix[int(cell)].split(',')
            all_cells_D.append(D_cellgen)
        tp_D[tp] = all_cells_D

    for tp, cell_gs in tp_D.items():
        noOfCells = len(cell_gs)
        noOfMut = len(cell_gs[0])
        #print(" Timepoint ",tp," cells ",noOfCells," Mutation ",noOfMut)
        missingEntries = 0
        for cg in cell_gs:
            for i in range(len(cg)):
                if cg[i] == '3':
                    missingEntries = missingEntries + 1
        #print(" Total missing entries ",missingEntries)
        missingRate = missingEntries / (noOfCells * noOfMut)
        print(" Timepoint ",tp," Missing rate ",missingRate)
        tp_MR[tp] = missingRate

    return tp_MR

''' First get the voted consensus genotype with unknowns. '''
def votedCG_WithUnknowns(cluster_cells_genotype):
    # Before you start this use the above function to calculate the missing rate
    cell_cg_len = len(cluster_cells_genotype)
    cg_len = len(cluster_cells_genotype[0])
    #print(" Inside Voted cluster genotype with unknowns ")
    #print(" Number of cluster cells ",cell_cg_len," mutations ",cg_len)

    #print(cell_gs_list[0])
    #print("============================")
    voted_cg = []
    #total0s = 0
    #total1s = 0
    # j = mutations, i = cells
    for j in range(0,cg_len):
        total0s = 0
        total1s = 0
        total3s = 0
        #print(" Vote cluster mutation ",j)
        for i in range(0,cell_cg_len):
            if cluster_cells_genotype[i][j] == '3':
                total3s = total3s+1
            if cluster_cells_genotype[i][j] == '0':
                total0s = total0s+1
            if cluster_cells_genotype[i][j] == '1':
                total1s = total1s+1
        #print(j," Total 0s ",total0s," Total 1s ",total1s," Total 3s ",total3s)
        #if total3s > total0s or total3s > total1s or total3s == total1s or total3s == total0s or total0s == total1s:
        #    voted_cg.append("unk")
        if total0s > (total1s + total3s): # Deterministic approach
            voted_cg.append('0')
        elif total1s > (total0s + total3s):
            voted_cg.append('1')
        else:
            voted_cg.append("unk")

    return voted_cg

''' Returns the value of Pr(c_1j^k, c_2j^k, …, _nj^k | cg_kj) based on the value of cell and cluster consensus genotype. '''
def prob_cell_j_k(cellMut, clusterMut, alpha, beta, missingRate):
    if (cellMut == '3' and clusterMut == "unk") or (cellMut == '1' and clusterMut == "unk") or (cellMut == '0' and clusterMut == "unk"):
        return 0.33
    if (cellMut == '3' and clusterMut == '1') or (cellMut == '3' and clusterMut == '0'):
        return missingRate
    if cellMut == '1' and clusterMut == '1':
        prob = (1 - missingRate) * (1 - beta)
        return prob
    if cellMut == '0' and clusterMut == '1':
        prob = (1 - missingRate) * beta
        return prob
    if cellMut == '1' and clusterMut == '0':
        prob = (1 - missingRate) * alpha
        return prob
    if cellMut == '0' and clusterMut == '0':
        prob = (1 - missingRate) * (1 - alpha)
        return prob

''' Get the prior value based on value of consensus genotype of subclone k for mutation j. '''
def prior_cg_kj(clusterMut, missingRate):
    if clusterMut == "unk":
        return missingRate
    if clusterMut == '1' or clusterMut == '0':
        return (1 - missingRate) / 2

''' Get the normalized value or P(D_j). The input will be cells likelihood cells_cj_cgkj, 
cluster mutation value, probability of all cells mutation given cluster mutation cgkj_cij 
and missingRate to calculate the prior. '''
def normalized_P_Dj(cells_cj_cgkj, clusterMut, cgkj_cij, missingRate):
    # For each mutation we will have the prob_cell_j_k and prior_cg_kj.
    # We just sum up all possible combinations when cg = 0, 1 or unk.
    P_Dj = 0
    possible_cg_kj = ['0','1','unk']

    for cg_kj in possible_cg_kj:
        if cg_kj == clusterMut: # we use this to not re-calculate this value
            P_cgkj_cij = cgkj_cij
        else:
            P_cgkj_cij =  cells_cj_cgkj * prior_cg_kj(cg_kj, missingRate)
        P_Dj = P_Dj + P_cgkj_cij

    return P_Dj

''' Calculate the probability of the consensus genotypes. '''
def probOf_votedCG(cluster_cells_tp, D_matrix, tp_alpha, tp_beta, tp_missingRate):
    cluster_prob = {}
    # First get the votedCG_WithUnknowns to get the voted CG for each subclones.
    voted_CG_unknowns = {}
    for cluster, cells in cluster_cells_tp.items():
        tp = cluster.split('_')[1]
        #print(" For cluster ",cluster," no. of cells ",len(cells))
        all_cells_D = []
        for cell in cells:
            D_cellgen = D_matrix[cell].split(',')
            all_cells_D.append(D_cellgen)
        
        #print(" All cells genotype D ",all_cells_D)
        voted_cluster_genotype = votedCG_WithUnknowns(all_cells_D)
        voted_CG_unknowns[cluster] = voted_cluster_genotype
        #print(" Voted cluster genotype ",voted_cluster_genotype)

    # Once you have the voted CG calculate Pr(cg_k | c_k) for each subclone with its supporting cells
    # Make sure to select appropriate alpha, beta and MR for each timepoint.
    for cluster, cells in cluster_cells_tp.items():
        tp = cluster.split('_')[1]
        voted_CG = voted_CG_unknowns[cluster]
        total_mutations = len(voted_CG)
        alpha = tp_alpha[tp]
        beta = tp_beta[tp]
        missingRate = tp_missingRate[tp]
        print(" For cluster ",cluster," voted cg ",voted_CG)
        sum_log_P_cgk_ck = 0
        # j = mutations, i = cell, k = subclones
        for j in range(len(voted_CG)):
            cells_cj_cgkj = 1 # cj = all cells of jth mutation, cgkj = consensus genotype of jth mutation of subclone k 
            for i in cells:
                D_cellgen = D_matrix[i].split(',')
                # cij = each individual cell i
                P_cij_cgkj = prob_cell_j_k(D_cellgen[j], voted_CG[j], alpha, beta, missingRate)
                cells_cj_cgkj = cells_cj_cgkj * P_cij_cgkj # multiply the likelihood of the cells

            prior_cgkj = prior_cg_kj(voted_CG[j], missingRate)
            #print(" P_cij_cgkj of all cells ",cells_cj_cgkj)
            #print(" Prior cgkj ",prior_cgkj)

            Pcgkj_cj = cells_cj_cgkj * prior_cgkj
            #print(" Numerator ",Pcgkj_cj)

            P_Dj = normalized_P_Dj(cells_cj_cgkj, voted_CG[j], Pcgkj_cj, missingRate)
            #print(" Normalizer P_Dj ",P_Dj)

            if P_Dj == 0:
                P_cgk_ck = 0
            else:
                P_cgk_ck = Pcgkj_cj / P_Dj
            #print(j," mutation prob ",P_cgk_ck)
            sum_log_P_cgk_ck = sum_log_P_cgk_ck + np.log(P_cgk_ck)
        
        #print(" Sum prob ",sum_log_P_cgk_ck)
        avg_sum_log_P_cgk_ck = sum_log_P_cgk_ck / total_mutations
        print(cluster," Cluster probability ",avg_sum_log_P_cgk_ck," # of cells ",len(cells))
        print(" =============================================================== ")
        cluster_prob[cluster] = avg_sum_log_P_cgk_ck
    return cluster_prob

''' Update the subclones to eliminate those below a certain threshold. '''
def updateSubclones_withThreshold(drop_cluster, cluster_cells_tp, tp_unassigned_cells):
    # With every cluster drop the unassigned cells will also change
    tp = drop_cluster.split('_')[1]
    if tp in tp_unassigned_cells:
        temp_cells = tp_unassigned_cells[tp]
        temp_cells.extend(cluster_cells_tp[drop_cluster])
        tp_unassigned_cells[tp] = temp_cells
    else:
        tp_unassigned_cells[tp] = cluster_cells_tp[drop_cluster]

    # After adding the unassigned cells we remove that cluster from cluster_cells_tp
    del(cluster_cells_tp[drop_cluster])

    print(" Unassigned cells ",tp_unassigned_cells)
    return cluster_cells_tp, tp_unassigned_cells

''' Select only those subclones < upper bound and sort them increasingly with probabilities. '''
def sortClusterProb(upper_bound, cluster_prob):
    selected_clusters = {}
    for cluster, prob in cluster_prob.items():
        if prob < upper_bound:
            selected_clusters[cluster] = prob

    # Then sort the selected clusters with increasing probability values
    sort_selected_clusters = sorted(selected_clusters.items(), key=lambda x:x[1])
    sort_selected_clusters = dict(sort_selected_clusters)
    sort_selected_clusters = list(sort_selected_clusters.keys())
    print(" Cluster prob ",cluster_prob)
    print(" Sort selected clusters ",sort_selected_clusters)
    return sort_selected_clusters

''' Get weighted thresholds for all timepoints. '''
def weighted_thresholds(cluster_prob, cluster_cells_tp):
    tp_thresholds = {}
    tp_cells = {} # Get the no. of cells at each timepoint
    tp_cluster_prob = {}

    for cluster, prob in cluster_prob.items():
        tp = cluster.split('_')[1]
        cells = cluster_cells_tp[cluster]
        cluster_prob = prob * len(cells)

        if tp in tp_cluster_prob:
            temp_cluster_prob = tp_cluster_prob[tp]
            tp_cluster_prob[tp] = temp_cluster_prob + cluster_prob
        else:
            tp_cluster_prob[tp] = cluster_prob

        if tp in tp_cells:
            temp_cells = tp_cells[tp]
            tp_cells[tp] = temp_cells + len(cells)
        else:
            tp_cells[tp] = len(cells)
    
    for tp, noOfCells in tp_cells.items():
        tp_prob = tp_cluster_prob[tp] / noOfCells
        #tp_thresholds[tp] = round(tp_prob,1)
        tp_thresholds[tp] = tp_prob

    print(" Timepoint prob ",tp_cluster_prob)
    print(" Timepoint cells ",tp_cells)
    print(" Timepoint thresholds ",tp_thresholds)
    # Find the max weighted avg of all time points which will be the upper bound. Return the upper bound.
    upper_bound = max(tp_thresholds.values())
    print(" Upper bound of the threshold ",upper_bound)
    return upper_bound

''' This is to maintain unassigned cells in each timepoint. '''
def assign_Empty_timepointCells(tp_missingRate):
    tp_prev_unassigned_cells = {}
    for tp in tp_missingRate:
        tp_prev_unassigned_cells[str(tp)] = []
    return tp_prev_unassigned_cells

''' Get parallel mutation count. '''
def getParallelMutCount(new_mutations):
    # 1, ..., m, for each mutation, count how many edges that have it as a new mutation not in the parent node. 
    # Minus that number by 1, and add up all the numbers for all mutations.
    mutation_node = {} # Get a dictionary with mutation as key and node as values
    for node, mutations in new_mutations.items():
        for mut in mutations:
            if mut in mutation_node:
                node_list = mutation_node[mut]
                node_list.append(node)
                mutation_node[mut] = node_list
            else:
                mutation_node[mut] = [node]
    print(" Mutations dict ",mutation_node)
    total_parallelMut = 0
    for mut, nodes in mutation_node.items():
        total_parallelMut = total_parallelMut + (len(nodes) - 1)

    print(" Total parallel mutation ",total_parallelMut)
    return total_parallelMut

''' Get back mutation count. '''
def getBackMutCount(all_mutations, new_mutations, parent_child):
    # back mutations are ones where a mutation present in parent is absent in child
    backMut = []
    print(" BACK MUTATIONS =================================== ")
    for parent, children in parent_child.items():
        if parent == -1:
            continue
        p_mut = all_mutations[parent]
        print(" Parent mut ",p_mut)
        for child in children:
            c_mut = all_mutations[child]
            print(" Child mut ",c_mut)
            c_new_mut = new_mutations[child]
            c_otherThan_newMut = set(c_mut) - set(c_new_mut)
            absent_mut = set(p_mut) - c_otherThan_newMut
            print(" Absent mut ",absent_mut)
            if absent_mut != set():
                backMut.extend(list(absent_mut))
    print(" Total back mutations ",len(backMut))
    print(" ================================================= ")
    return len(backMut)

''' Select the tree with least parallel and back mutations. '''
def selectBestTree(all_trees, cloneNodes_forEachRun, tp_cells_forEachRun):
    tree_parallel_back_count = {} # use this to save the count for each tree
    for i in range(len(all_trees)):
        # need a dict with nodeId as keys and new mutations for parallel mut count as values
        # need another dict with nodeId all mutations for back mutation. Another dict with parent as key and child as value.
        Graph = all_trees[i]
        new_mutations = {}
        all_mutations = {}
        parent_child = {}
        for k in range(len(Graph)):
            #up_Graph[i].id," ",up_Graph[i].timepoint_," ",up_Graph[i].parentID," ",up_Graph[i].type," ",up_Graph[i].mutations
            nodeId = Graph[k].id # this will also be a child node
            pId = (Graph[k].parentID)[0]
            mutations = Graph[k].mutations
            if pId in parent_child:
                temp_child = parent_child[pId]
                temp_child.append(nodeId)
                parent_child[pId] = temp_child
            else:
                parent_child[pId] = [nodeId]

            all_mutations[nodeId] = mutations

        #print(" All mutations ",all_mutations)
        #print(" Parent child ",parent_child)
        # Next we get the new_mutations by comparing parent and child mutations
        for parent, child in parent_child.items():
            if parent == -1:
                continue
            p_mut = all_mutations[parent]
            for c in child:
                c_mut = all_mutations[c]
                new_mutations[c] = set(c_mut) - set(p_mut)

        print(" New mutations ",new_mutations)
        total_parallelMut = getParallelMutCount(new_mutations)
        backMutations = getBackMutCount(all_mutations, new_mutations, parent_child)
        tree_parallel_back_count[i] = total_parallelMut + backMutations

    print(" All trees parallel and back mutation count ",tree_parallel_back_count)
    treeIndex_ForMin_mutations = min(tree_parallel_back_count, key=tree_parallel_back_count.get)
    print(" Tree index with least no. of parallel and back mutations ",treeIndex_ForMin_mutations)
    return all_trees[treeIndex_ForMin_mutations], cloneNodes_forEachRun[treeIndex_ForMin_mutations], tp_cells_forEachRun[treeIndex_ForMin_mutations]

''' Calculate cell fraction for each subclones. '''
def calculate_cellFraction(tp_cell, original_cluster_cells_tp):
    cluster_cellFraction = {}
    for cluster, cells in original_cluster_cells_tp.items():
        tp = cluster.split('_')[1]
        total_cells = len(tp_cell[tp])
        cluster_cellFraction[cluster] = len(cells)/total_cells
    print(" Cluster cell fraction ",cluster_cellFraction)
    return cluster_cellFraction

''' Update search space to include only subclones with <= 3 cells and cell fraction < 0.1. '''
def updateSearchSpace(clusters_to_eliminate, original_cluster_cells_tp, cluster_cellFraction):
    updated_clusters_to_eliminate = []
    for small_clones in clusters_to_eliminate:
        originalCells = len(original_cluster_cells_tp[small_clones]) # extra checks before eliminating clones
        clusterCellFraction = cluster_cellFraction[small_clones]
        print(" Clone ",small_clones," originalCells ",originalCells," cell fraction ",clusterCellFraction)
        if originalCells <= 3 and clusterCellFraction < 0.1:
            updated_clusters_to_eliminate.append(small_clones)
    print(" Updated search space ",updated_clusters_to_eliminate)
    return updated_clusters_to_eliminate

''' Get the tree with a search space where some clones needs to be eliminated. '''
def eliminateClonesAndInferTree(clusters_to_eliminate, cluster_cells_tp, tp_prev_unassigned_cells, cells_genotype, cluster_genotype):
    graph_forEachRun = [] # Save the graph after eliminating all subclones
    cloneNodes_forEachRun = [] # Save the clone nodes after each run
    tp_cells_forEachRun = [] # Save the no. of cells assigned to each cluster
    for small_clones in clusters_to_eliminate:
        print("========================= NEW ITERATION ================================ ")
        cluster_cells_tp, tp_prev_unassigned_cells = updateSubclones_withThreshold(small_clones, cluster_cells_tp, tp_prev_unassigned_cells)
        print(" After dropping cluster ",small_clones)
        print(" Updated clusters in timepoints ",cluster_cells_tp.keys())
        print(" TP unassigned cells ",tp_prev_unassigned_cells)

        print(" After recalculating FP FN ")
        tp_weighted_fp, tp_weighted_fn, voted_D_cg, cells_ignore_mutations = getTimepoint_FPFN(cluster_cells_tp, cells_genotype)
        print(" Updated clusters in timepoints ",cluster_cells_tp)

        # Below we use the tp_prev_unassigned_cells to make sure that these cells are not considered while getting the consensus genotype
        tp_clusters, tp_wrongClones, tp_true_clusters, tp_cell_cluster_prob, tp_unassigned_cells = get_WrongClones(cells_genotype, voted_D_cg, tp_cell, cluster_cells_tp, cells_ignore_mutations, tp_weighted_fp, tp_weighted_fn, tp_prev_unassigned_cells)
        print(" After assignment of cells in clusters ")
        print(" After wrong clones unassigned cells ")

        for tp, cells in tp_unassigned_cells.items():
            print(" Timepoint ",tp," no. of cells ",len(cells))

        print(" True cluster assignment is empty ",tp_true_clusters)
        print(" After dropping cluster ",small_clones)
        print(" Updated clusters in timepoints ",cluster_cells_tp.keys())
        print(" Cell cluster prob ",tp_cell_cluster_prob[3].keys())

        for tp, cluster_cell in tp_true_clusters.items():
            print(" Timepoint ",tp)
            for cluster, cell in cluster_cell.items():
                print(" Cluster ",cluster," # of cells ",len(cell))

        correced_tp_clusters, tp_true_clusters = updateSubclones(tp_true_clusters, tp_unassigned_cells, tp_prev_unassigned_cells, tp_cell_cluster_prob, cells_genotype)
        #if correced_tp_clusters == "Don't check further": # Because all cells are unassigned and no cluster is selected in this timepoint
        #    continue
        print(" Corrected cluster genotype ",correced_tp_clusters)

        # Merge those clusters having same mutations or consensus genotype.
        correced_tp_clusters, tp_true_clusters = merge_tp_clusters(correced_tp_clusters, tp_true_clusters)
        #correced_tp_clusters, tp_true_clusters = merge_clones_WithOneCell(correced_tp_clusters, tp_true_clusters)
        tp_cluster_cells = {}
        print(" After reassigning and merging cells ")
        for tp, cluster_cell in tp_true_clusters.items():
            print(" Timepoint ",tp)
            for cluster, cell in cluster_cell.items():
                print(" Cluster ",cluster," # of cells ",len(cell))
                tp_cluster_cells[cluster+"_"+str(tp)] = len(cell)

        print(" Tp cluster cell dict ",tp_cluster_cells)
        print(" Inital Graph with subclones placed in each timepoint ")
        Graph, node_tp, clone_node = build_graph(correced_tp_clusters)
        print(" Clone to nodes ",clone_node)
        #mergeClones(tp_clusters, tp_wrongClones, cluster_cells_tp, cells_genotype, cluster_genotype, args.alpha, args.beta)
        # =============================================
        # Correcting clusters ends here.
        # Longitudinal algorithm to infer the tree after clustering starts here.
        # We first get one layer of unobserved subclones
        up_Graph, usc_children = find_unobservedSubclones(Graph, node_tp, False) # Updates the graph by connecting subclones with same mutations in each timepoint, adding edges and branch length.
        print(" Timepoint usc ",usc_children)
        up_Graph, node_tp = select_unobservedSubclones(up_Graph, node_tp, usc_children, False)
        print(" Check for second level of unobserved subclones ")
        print(" ============================================== ")
        up_Graph, second_usc_children = find_unobservedSubclones(up_Graph, node_tp, True)
        print(" Timepoint second usc ",second_usc_children)
        up_Graph, node_tp = select_unobservedSubclones(up_Graph, node_tp, second_usc_children, True)

        # Traverse the updated graph with unobserved subclones and find unconnected nodes.
        up_Graph = find_nodeV(up_Graph, node_tp)
        graph_forEachRun.append(up_Graph)

        for i in range(len(up_Graph)):
            print(up_Graph[i].id," ",up_Graph[i].timepoint_," ",up_Graph[i].parentID," ",up_Graph[i].type," ",up_Graph[i].mutations," ",up_Graph[i].edge_length)

        cloneNodes_forEachRun.append(clone_node)
        tp_cells_forEachRun.append(tp_cluster_cells)
    return graph_forEachRun, cloneNodes_forEachRun, tp_cells_forEachRun

        
parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input Gprime matrix (Input to the Longitudinal tree)")
parser.add_argument("-cg", "--cg",dest ="cg", help="Input GDoublePrime matrix")
parser.add_argument("-D", "--D",dest ="D", help="D matrix representing true cell genotypes.")
parser.add_argument("-cc", "--cc",dest ="cc", help="Assignment file indicating cells belonging to clusters.")
parser.add_argument("-celltp", "--celltp",dest = "celltp", help="Cell timepoints.")
#parser.add_argument("-alpha","--alpha",dest = "alpha", help="False positive rate")
#parser.add_argument("-beta","--beta",dest = "beta", help="False negative rate")
#parser.add_argument("-FPFN","--FPFN",dest = "FPFN", help="Estimated False positive and False negative rate from BnpC.")
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
print(" Cluster cell ")
print(cluster_cell)
print(" ============ ")
print(" Cluster mutation ")
cluster_mutation = get_cluster_mutations(args.input)
print(cluster_mutation)
print(" ============ ")
cell_tp, tp_cell = get_cell_timepoints(args.celltp)
cluster_cells_tp = getCellsOfEachCluster(cell_tp, cluster_cell)
print(" Cell timepoints ",tp_cell)
print(" Original cluster cells timepoint ")
print(cluster_cells_tp)
original_cluster_cells_tp = cluster_cells_tp # use this for search space refining
# Now from original_cluster_cells_tp and tp_cell calculate cell fraction
cluster_cellFraction = calculate_cellFraction(tp_cell, original_cluster_cells_tp)

# Correcting clusters starts here. 
# ============================================
cells_genotype = getGenotype_dict(args.D, "D") # DMatrix
cells_G_genotype = getGenotype_dict(args.cg, "G") # Get the cell genotype from the consensus genotype inferred by BnpC.
#print(" D matrix ")
#print(cells_genotype)
#print(" Cells genotype ")
#print(cells_G_genotype)
tp_weighted_fp, tp_weighted_fn, voted_D_cg, cells_ignore_mutations = getTimepoint_FPFN(cluster_cells_tp, cells_genotype)
tp_missingRate = timepoint_missingRate(tp_cell, cells_genotype)
# we calculate the fidelity (existence) of each consensus genotype
cluster_prob = probOf_votedCG(cluster_cells_tp, cells_genotype, tp_weighted_fp, tp_weighted_fn, tp_missingRate)

# Get weighted thresholds for each timepoint
upper_bound = weighted_thresholds(cluster_prob, cluster_cells_tp)

# The following clusters are to be checked for elimination
clusters_to_eliminate = sortClusterProb(upper_bound, cluster_prob)
tp_prev_unassigned_cells = assign_Empty_timepointCells(tp_missingRate) # This is just to maintain unassigned cells in timepoints

cluster_genotype = getGenotype_dict(args.input, "cluster")
print(" BnpC cluster genotype ")
print(cluster_genotype)

graph_forEachRun = [] # Save the graph after eliminating all subclones
cloneNodes_forEachRun = [] # Save the clone nodes after each run
tp_cells_forEachRun = [] # Save the no. of cells assigned to each cluster

# update clusters_to_eliminate to only include subclones with no. of cells <= 3 and cell fraction < 0.1
clusters_to_eliminate = updateSearchSpace(clusters_to_eliminate, original_cluster_cells_tp, cluster_cellFraction)
if clusters_to_eliminate != []:
    # Loop over the search space and eliminate the subclones
    graph_forEachRun, cloneNodes_forEachRun, tp_cells_forEachRun = eliminateClonesAndInferTree(clusters_to_eliminate, cluster_cells_tp, tp_prev_unassigned_cells, cells_genotype, cluster_genotype)
else:
    # Since there are no subclones to eliminate we get only one tree 
    # Below we use the tp_prev_unassigned_cells to make sure that these cells are not considered while getting the consensus genotype 
    tp_clusters, tp_wrongClones, tp_true_clusters, tp_cell_cluster_prob, tp_unassigned_cells = get_WrongClones(cells_genotype, voted_D_cg, tp_cell, cluster_cells_tp, cells_ignore_mutations, tp_weighted_fp, tp_weighted_fn, tp_prev_unassigned_cells)
    print(" After assignment of cells in clusters ")
    print(" After wrong clones unassigned cells ")

    for tp, cells in tp_unassigned_cells.items():
        print(" Timepoint ",tp," no. of cells ",len(cells))
    #print(" Clusters to cell dict ",cluster_cell)
    #print(" Wrong clones ",tp_wrongClones)
    #print(" True clusters ",tp_true_clusters)
    #print(" Timepoint clusters ",tp_clusters)
    #print(" Cell cluster prob ",tp_cell_cluster_prob)
    print(" Before reassigning cells ")
    #print(" TP true clusters ",tp_true_clusters)

    #if tp_true_clusters == {}:
    print(" True cluster assignment is empty ",tp_true_clusters)
    #print(" After dropping cluster ",small_clones)
    print(" Updated clusters in timepoints ",cluster_cells_tp.keys())
    print(" Cell cluster prob ",tp_cell_cluster_prob[3].keys())

    for tp, cluster_cell in tp_true_clusters.items():
        print(" Timepoint ",tp)
        for cluster, cell in cluster_cell.items():
            print(" Cluster ",cluster," # of cells ",len(cell))

    
    correced_tp_clusters, tp_true_clusters = updateSubclones(tp_true_clusters, tp_unassigned_cells, tp_prev_unassigned_cells, tp_cell_cluster_prob, cells_genotype)
    #if correced_tp_clusters == "Don't check further": # Because all cells are unassigned and no cluster is selected in this timepoint
    #    continue
    print(" Corrected cluster genotype ",correced_tp_clusters)

    # Merge those clusters having same mutations or consensus genotype.
    correced_tp_clusters, tp_true_clusters = merge_tp_clusters(correced_tp_clusters, tp_true_clusters)

    #correced_tp_clusters, tp_true_clusters = merge_clones_WithOneCell(correced_tp_clusters, tp_true_clusters)
    tp_cluster_cells = {}
    print(" After reassigning and merging cells ")
    for tp, cluster_cell in tp_true_clusters.items():
        print(" Timepoint ",tp)
        for cluster, cell in cluster_cell.items():
            print(" Cluster ",cluster," # of cells ",len(cell))
            tp_cluster_cells[cluster+"_"+str(tp)] = len(cell)

    print(" Tp cluster cell dict ",tp_cluster_cells)
    print(" Inital Graph with subclones placed in each timepoint ")
    Graph, node_tp, clone_node = build_graph(correced_tp_clusters)
    print(" Clone to nodes ",clone_node)
    #mergeClones(tp_clusters, tp_wrongClones, cluster_cells_tp, cells_genotype, cluster_genotype, args.alpha, args.beta)
    # =============================================
    # Correcting clusters ends here.

    # Longitudinal algorithm to infer the tree after clustering starts here.
    # We first get one layer of unobserved subclones 
    up_Graph, usc_children = find_unobservedSubclones(Graph, node_tp, False) # Updates the graph by connecting subclones with same mutations in each timepoint, adding edges and branch length.
    print(" Timepoint usc ",usc_children)
    up_Graph, node_tp = select_unobservedSubclones(up_Graph, node_tp, usc_children, False)
    print(" Check for second level of unobserved subclones ")
    print(" ============================================== ")
    up_Graph, second_usc_children = find_unobservedSubclones(up_Graph, node_tp, True)
    print(" Timepoint second usc ",second_usc_children)
    up_Graph, node_tp = select_unobservedSubclones(up_Graph, node_tp, second_usc_children, True)

    # Traverse the updated graph with unobserved subclones and find unconnected nodes.
    up_Graph = find_nodeV(up_Graph, node_tp)
    graph_forEachRun.append(up_Graph)

    for i in range(len(up_Graph)):
        print(up_Graph[i].id," ",up_Graph[i].timepoint_," ",up_Graph[i].parentID," ",up_Graph[i].type," ",up_Graph[i].mutations," ",up_Graph[i].edge_length)

    cloneNodes_forEachRun.append(clone_node)
    tp_cells_forEachRun.append(tp_cluster_cells)

# All the graphs after the iteration
print(" All graphs from each iteration ",len(graph_forEachRun))
bestGraph, bestCloneNode, bestClusterCells = selectBestTree(graph_forEachRun, cloneNodes_forEachRun, tp_cells_forEachRun)
print(" BEST Graph ========================================== ")
for i in range(len(bestGraph)):
    print(bestGraph[i].id," ",bestGraph[i].timepoint_," ",bestGraph[i].parentID," ",bestGraph[i].type," ",bestGraph[i].mutations," ",bestGraph[i].edge_length)

# We will save the final graph here. Need to modify the save_graph function for that.
#save_graph(up_Graph,"default_mr_op/rep1_tree.csv")
save_graph(bestGraph,args.op_tree)
np.save(args.cloneNode, bestCloneNode) # Save the clone to nodes mapping to be used for evaluation
np.save(args.cloneCells, bestClusterCells) # Save the no of cells in each clones
