#!/usr/bin/env python2
# -*- coding: utf-8 -*-
"""
@authors: Xian Fan Mallory
Contacting email: fan@cs.fsu.edu
"""

import sys
import argparse
import numpy as np
import random
import copy

class Edge():
    def __init__(self, p, c):
        self.p = p
        self.c = c

class Node():
    def __init__(self, e, p, c):
        # only the edge above this node is listed here
        self.e = e
        self.p = p
        # this node may have multiple children, c is an array
        self.c = c

# initialize a matrix that has r rows and c columns
def init_m(r, c, value):
    m = []
    for i in range(r):
        m.append(copy.deepcopy([value]*c))
    return m



# return a n*m matrix with zero entries
def init(n, m):
    ret = []
    for i in range(n):
        ret.append([0]*m)
        for j in range(m):
            ret[i][j] = 0
    return ret

   
def print_matrix(M, n, matrix_file):
    matrix_f = open(matrix_file, "w")
    for i in range(n):
        str_ = [str(j) for j in M[i]]
        print("\t".join(str_), file = matrix_f) 
    matrix_f.close()

def print_f(array, out_file):
    out_f = open(out_file, 'w')
    for i in range(len(array)):
        print("\t".join(array[i]), file = out_f)
    out_f.close()

# return an array of mutations that are above an edge
def mut_above_edge(edgeID, e_dict, n_dict, mut_array):
    # turn mut_array into a dict. mut_array has edgeID and mut joined by semicolon
    mut_dict = {}
    for i in range(len(mut_array)):
        mut_dict[mut_array[i][0]] = mut_array[i][1]
    p = e_dict[edgeID].p
    ret = []
    if mut_dict[edgeID] != "NA":
        ret = mut_dict[edgeID].split(";")
    while p != "-1":
        edgeID = n_dict[p].e
        if mut_dict[edgeID] != "NA":
            for i in mut_dict[edgeID].split(";"):
                ret.append(i)
        p = e_dict[edgeID].p
    #print(str(ret))
    return ret

# return the array of leaf IDs under this edge
def leaves_under_edge(edgeID, e_dict, n_dict):
    ret_node = []
    node = e_dict[edgeID].c
    # a stack to store all the nodes below for retrieval of the leaves
    c_arr = [node]
    while c_arr != []:
        node = c_arr.pop()
        #print(node)
        if n_dict[node].c != "NA":
            # more children
            for i in n_dict[node].c.split(";"):
                c_arr.append(i)
        else:
            ret_node.append(node)

    return ret_node

# given a leaf ID, and an edge dict that has edge ID as the key, .c and .p as the child and parent node ID, return all the edges above this leaf ID in array 
def retrieve_edges(leafID, n_dict):
    p = n_dict[leafID].p
    e = n_dict[leafID].e
    e_array = [e]
    while p != "-1":
        ID = p
        p = n_dict[ID].p
        e = n_dict[ID].e
        e_array.append(e)
    return e_array


def read_edge_num(tree_f):
    ret = 0
    file = open(tree_f, "r")
    line = file.readline()
    while(line != ""):
        ret = ret + 1
        line = file.readline()
    file.close()
    return ret 

# If no mutation on this edge, then "NA" on the second column
def distribute_mutations(mut_n, tree_f, out_f):
    # to make all mutations add up to mut_n
    edge_num = read_edge_num(tree_f)
    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    mut_array = []
    # mutation ID starts from 0 (0-based)
    prev = -1
 
    # to avoid negative number of mutation, go through the whole file once and if negative, go through each line again
    edge = []
    while(line != ""):
        line_a = line.split('\t')
        # all branch lengths have already been normalized. They add up to 1
        branch_len = float(line_a[3])
        edge.append([line_a[0], branch_len])
        line = file.readline().rstrip('\n')
    file.close()

    tag = True
    mut_nums = init_m(len(edge), 2, 0)
    while tag:
        prev = 0
        for it in range(len(edge)):
            edgeID, branch_len = edge[it]
            mut_num = round(mut_n * branch_len)

            if it == len(edge) - 1:
                mut_num = mut_n - prev
            mut_nums[it] = [edgeID, mut_num]
            if mut_num < 0:
                print("In distributing mutations, mut_num becomes negative. Start again. ")
                break
            prev = prev + mut_num
            if it == len(edge) - 1:
                print("Success in distributing mutations. ")
                tag = False

    prev = -1
    for i in range(len(mut_nums)):
        edgeID, mut_num = mut_nums[i]
        mut_IDs = "NA"
        if mut_num != 0:
            mut_IDs = str(prev + 1)
            for j in range(prev + 2, prev + mut_num + 1):
                mut_IDs = mut_IDs + ";" + str(j)
        mut_array.append([edgeID, mut_IDs])
        prev = prev + mut_num

    print_f(mut_array, out_f)
    return mut_array


def read_leaf_num(tree_f):
    ret = 0
    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    while(line != ""):
        line_a = line.split('\t')
        if line_a[5] == "1":
            ret = ret + 1
        line = file.readline().rstrip('\n')

    file.close()
    return ret 


def distribute_SNVcells(cell_n, tree_f, out_f, eta):
    # to make all cell number add up to cell_n
    leaf_num = read_leaf_num(tree_f)
    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    SNVcell_array = []
    # SNV cell ID starts from 0 (0-based)
    prev = -1

    # to avoid negative number of cell, go through the whole file once and if negative, go through each line again
    leafnode = []
    while(line != ""):
        line_a = line.split('\t')
        # only look at leaves' percentages
        if line_a[5] == "1": 
            SNVcellP = float(line_a[4])
            leafnode.append([line_a[2], SNVcellP])
        line = file.readline().rstrip('\n')

    file.close()

    tag = True
    SNVcell_nums = init_m(len(leafnode), 2, 0)

    # simulate non doublet cells first
    while tag:
        prev = 0
        for it in range(len(leafnode)):
            leafID, SNVcellP = leafnode[it]
            SNVcell_num = round(cell_n * (1 - eta) * SNVcellP)
            if it == len(leafnode) - 1:
                SNVcell_num = cell_n * (1 - eta) - prev
            SNVcell_nums[it] = [leafID, SNVcell_num]
            if SNVcell_num < 0:
                print("In distributing SNV cells, SNVcell_num becomes negative. Start again. ")
                break
            prev = prev + SNVcell_num
            if it == len(leafnode) - 1:
                print("Success in distributing SNV cells. ")
                tag = False

    # simulate doublet cells
    # leafIDs (separated by ;) to cellIDs (separated by ;)
    dbl_dict = {}
    for i in range(cell_n - cell_n * (1 - eta)):
        # randomly pick two leaves
        leafIDs = sorted(random.sample(range(len(leafnode)), 2))
        leafIDs = [str(x) for x in leafIDs]
        cellID = str(i + cell_n * (1 - eta))
        if ";".join(leafIDs) not in dbl_dict.keys():
            dbl_dict[leafIDs] = cellID
        else:
            dbl_dict[leafIDs] = dbl_dict[leafIDs] + ";" + cellID 

    prev = -1
    for i in range(len(SNVcell_nums)):
        leafID, SNVcell_num = SNVcell_nums[i]
        SNVcell_IDs = "NA"
        if SNVcell_num != 0:
            SNVcell_IDs = str(prev + 1)
            for j in range(prev + 2, prev + SNVcell_num + 1):
                SNVcell_IDs = SNVcell_IDs + ";" + str(j)
        SNVcell_array.append([leafID, SNVcell_IDs])
        prev = prev + SNVcell_num

    # add dbl_dict to SNVcell_array
    for i in dbl_dict.keys():
        SNVcell_array.append([i, dbl_dict[i]])

    print_f(SNVcell_array, out_f)
    return SNVcell_array


def add_missing(G, n, m, missingP):
    # in G (n*m) matrix, select missingP * n * m entries and flip them to 3. 
    missing_arr = random.sample(range(n * m), int(n * m * missingP))
    for i in missing_arr:
        row = int(i / m)
        column = i % m
        G[row][column] = 3
    return G


def count_total_value(M, n, m, value):
    total = 0
    for i in range(n):
        for j in range(m):
            if M[i][j] == value:
                total = total + 1
    return total


def add_FPFNs(D, n, m, alpha, beta, total_zeros, total_ones):
    FP_arr = random.sample(range(total_zeros), int(total_zeros * alpha))
    FN_arr = random.sample(range(total_ones), int(total_ones * beta))
    D_ = copy.deepcopy(D)
    index_neg = 0
    index_pos = 0
    for i in range(n):
        for j in range(m):
            if D[i][j] == 0:
                if index_neg in FP_arr:
                    # flip it
                    D_[i][j] = 1
                index_neg = index_neg + 1
            elif D[i][j] == 1:
                if index_pos in FN_arr:
                    # flip it
                    D_[i][j] = 0
                index_pos = index_pos + 1
    return D_
    

# from tree_f we know which cells (on leaf nodes) have which mutations (on edges)
def mutation_matrix(mut_array, SNVcell_array, tree_f, n, m, missingP, alpha, beta, G_matrix_file, D_matrix_file, D_miss_file):

    # read mut_array to a dict whose first column is the edge IDs and second are the mutation IDs 
    mut_array_dict = {}
    for i in range(len(mut_array)):
        mut_array_dict[mut_array[i][0]] = mut_array[i][1]

    file = open(tree_f, "r")
    line = file.readline().rstrip('\n')
    # for each leaf node, record all the edges that are above it
    # for each edge, record all the cells are below it


    # node_dict: a dict that has the child node ID as the key, .e and .p as the values, the edge above it and the parent node ID, respectively
    # edge_dict: a dict that has the edge ID as the key, .c and .p as the values, the child and parent node IDs on the two ends, respectivley
    edge_dict = {}
    node_dict = {}
    # c is a temporary data structure for retrieving all children nodes for all parent nodes for node_dict
    c_dict = {}
    while(line != ""):
        line_a = line.split('\t')
        # line_a[1] is parent, line_a[2] is child
        edge_dict[line_a[0]] = Edge(line_a[1], line_a[2])
        node_dict[line_a[2]] = Node(line_a[0], line_a[1], "NA")
        # record the children for each parent to be added to node_dict
        if line_a[1] not in c_dict.keys():
            c_dict[line_a[1]] = line_a[2]
        else:
            c_dict[line_a[1]] = c_dict[line_a[1]] + ";" + line_a[2]
        line = file.readline().rstrip('\n')

    file.close()

    for i in node_dict.keys():
        if i in c_dict.keys():
            node_dict[i].c = c_dict[i]
        else:
            # this is a leaf
            node_dict[i].c = "NA"
 
    # a dict from SNV cell ID to the mutation ID array
    SNVcell_mut = {}

    # cell to mutation ID arrays
    for i in range(len(SNVcell_array)):
        leafID = SNVcell_array[i][0]
        SNVcellIDs = SNVcell_array[i][1]
        # it's possible no SNV cells is assigned to this leaf
        if SNVcellIDs == "NA":
            continue
        # given the leaf ID, return the edges above this leaf in array
        # for doublets, consider all edges for both paths
        leafIDs = leafID.split(";")
        edges_above_leaf = []
        for l in leafIDs:
            edges_above_leaf.append(retrieve_edges(l, node_dict))

        # mutation IDs in an array corresponding to this SNV cell
        mutID_array = []
        for e in edges_above_leaf:
            if mut_array_dict[e] != "NA":
                for mutID in mut_array_dict[e].split(";"):
                    mutID_array.append(mutID)
   
        for j in SNVcellIDs.split(";"):
            SNVcell_mut[j] = copy.deepcopy(mutID_array)

    # construct G matrix. n*m: n(cell), m(mut)
    G = init(n, m)
    for cellID in range(n):
        if str(cellID) in SNVcell_mut.keys():
            mut_array_ = SNVcell_mut[str(cellID)]
            for mutID in range(m):
                if str(mutID) in mut_array_:
                    G[cellID][mutID] = 1

    print_matrix(G, n, G_matrix_file)
    total_zeros = count_total_value(G, n, m, 0)
    total_ones = count_total_value(G, n, m, 1)
    #print("total zeros for G: " + str(total_zeros) + "; total ones for G: " + str(total_ones))

    # Step a. now make missing data to produce D 
    D_miss = add_missing(G, n, m, missingP) 
    print_matrix(D_miss, n, D_miss_file)

    # Step b, c. flip 0's by \alpha and 1's by \beta 
    total_zeros = count_total_value(D_miss, n, m, 0)
    total_ones = count_total_value(D_miss, n, m, 1)
    #print("total zeros after missing: " + str(total_zeros) + "; total ones after missing: " + str(total_ones))
    D_miss_FP_FN = add_FPFNs(D_miss, n, m, alpha, beta, total_zeros, total_ones)
    total_zeros = count_total_value(D_miss_FP_FN, n, m, 0)
    total_ones = count_total_value(D_miss_FP_FN, n, m, 1)
    #print("total zeros after FPFN: " + str(total_zeros) + "; total ones after FPFN: " + str(total_ones))

    # store the D matrix
    print_matrix(D_miss_FP_FN, n, D_matrix_file)
    

   
# beginning of the program
if len(sys.argv) <= 1:
    print("""
    This generates the mutation matrix with the ground truth data. 
    Usage: python sim_par.py -a [alpha] -b [beta] -m [missing-rate] -c [num-cells] -n [num-mutations] -f [input-tree-file] -P [prefix-output-files]
        -a (--alpha)        False positive rate. [0.01]
        -b (--beta)         False negative rate. [0.2]
        -m (--missing-rate) Missing rate in G. [0.2]
        -c (--num-cells)    Total number of SNV cells. [500]
        -n (--num-mutations)    Total number of mutations. [200]
        -e (--doublet-rate) Doublet rate. [0]
        -f (--tree-file)    The input tree structure file. ["NA"]
        -P (--prefix)       Prefix of output files. 
    """)
    sys.exit(0)

parser = argparse.ArgumentParser(description='This script iterates all parameters to generate simulated data.')
parser.add_argument('-a', '--alpha', default=0.01)
parser.add_argument('-b', '--beta', default=0.2)
parser.add_argument('-m', '--missing-rate', default=0.2)
parser.add_argument('-c', '--num-cells', default=500)
parser.add_argument('-n', '--num-mutations', default=200)
parser.add_argument('-e', '--doublet-rate', default=0)
parser.add_argument('-f', '--tree-file', default="NA")
parser.add_argument('-P', '--prefix', default="NA")


args = parser.parse_args()
alpha = float(args.alpha)
beta = float(args.beta)
missing = float(args.missing_rate)
cell_n = int(args.num_cells)
mut_n = int(args.num_mutations)
eta = float(args.doublet_rate)
tree_f = args.tree_file
prefix = args.prefix

# First, distribute the mutations on the edges. Read tree_f output to prefix + ".mut.csv" 
mut_array = distribute_mutations(mut_n, tree_f, prefix + ".mut.csv")

# Second, distribute the cells on the leaf nodes. Read tree_f output to prefix + ".SNVcell.csv" 
SNVcell_array = distribute_SNVcells(cell_n, tree_f, prefix + ".SNVcell.csv", eta)

# Last, make mutation matrices G and D. 
mutation_matrix(mut_array, SNVcell_array, tree_f, cell_n, mut_n, missing, alpha, beta, prefix + ".G.csv", prefix + ".D.csv", prefix + ".miss.csv")


