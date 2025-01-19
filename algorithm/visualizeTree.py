# Author: Rituparna Khan

import pydot
import argparse
import random
import numpy as np

''' Get the list of mutations having genes. '''
# Input is list of mutations having genes and mutations present in the subclone
def mutGenes(mutList, nodeMut):
    mut = []
    for m in mutList:
        if m in nodeMut:
            mut.append(m)
    #mutGeneList = ";".join(str(i) for i in mut)
    return mut

''' Read the txt file to save small and large SA501 dict. '''
def getMutDict():
    filen = open('../large_smallMut.txt','r')
    fline = filen.readline().rstrip('\n')
    mutDict = {}
    while(fline != ''):
        key = fline.split("\t")[0]
        val = fline.split("\t")[1]
        mutDict[key] = val
        fline = filen.readline().rstrip('\n')
    #print("Mut dict ",mutDict)
    return mutDict

''' Map the small and large SA501 mut. '''
# Input is mutations of the current node
def getMutList(nodeMut):
    mutDict = getMutDict()
    smallMut = []
    for k,v in mutDict.items():
        if int(k) in nodeMut:
            smallMut.append(v)
    #print("Small Mut ",smallMut)
    return smallMut

''' Use this method to map the gene mutations with its original no. mutations. '''
# Input is mutations of the current node
def gene_smallSA501(nodeMut):
    # This dict has gene mutations mapped to original mutation nos.
    gene_all = {0: 0, 1: 1, 2: 2, 3: 7, 4: 8, 5: 9, 6: 10, 7: 11, 8: 14, 9: 17, 10: 18, 11: 19, 12: 20, 13: 21, 14: 22, 15: 23, 16: 24, 17: 25, 18: 26, 19: 27, 20: 28, 21: 30, 22: 31, 23: 32, 24: 33, 25: 34, 26: 35, 27: 37, 28: 38, 29: 40, 30: 41, 31: 42, 32: 43, 33: 44, 34: 45, 35: 50, 36: 51, 37: 54}
    largeMut = []
    for m in nodeMut:
        largeMut.append(gene_all[m])
    return largeMut

def AML_mut(nodeMut):
    #AML38 = {0: 'NPM1_p.L287fs', 1: 'IDH1_p.R132H', 2: 'IDH2_p.R140Q', 3: 'FLT3-ITD', 4: 'FLT3_p.D835H', 5: 'PTPN11_p.D61H', 6: 'PTPN11_p.A72G', 7: 'PTPN11_p.G503A', 8: 'KRAS_p.G12A', 9: 'KRAS_p.G12D', 10: 'NRAS_p.G13R', 11: 'NRAS_p.G12A'}
    #AML99 = {0: 'RUNX1_p.R162K', 1: 'DNMT3A_p.R882H', 2: 'IDH2_p.R140Q', 3: 'CSF3R_p.Q768X', 4: 'SRSF2_p.P95L', 5: 'NRAS_p.G60E', 6: 'RUNX1_p.S318fs', 7: 'PTPN11_p.I282V', 8: 'FLT3-ITD', 9: 'IDH1_p.R132C', 10: 'CSF3R_p.Q776X'}
    AML09 = {0: 'NPM1_p.L287fs', 1: 'FLT3-ITD', 2: 'FLT3_p.D835E', 3: 'FLT3_p.D835Y', 4: 'KRAS_p.G13D', 5: 'WT1_p.P372fs'}
    mut = []
    for m in nodeMut:
        mut.append(AML09[m])
    return mut

def drawTree(tp_nodes, id_cells, id_newMut, id_plMut, id_backMut, id_edges, usc_nodes, name, sample):
    print("Inside visualize Tree ")
    print("Back mutation ",id_backMut)
    # For small SA501 mutations having genes
    #mutList = [3,12,15,18,9,0,6,5,14,2,13,10,19,8,1,11]
    # For large SA501 mutations having genes
    #mutList = [0,1,2,7,8,9,10,11,14,17,18,19,20,21,22,23,24,25,26,27,28,30,31,32,33,34,35,37,38,40,41,42,43,44,45,50,51,54]

    # Node Ids should be cell nos. Edges should be mutations indicating both new and loss. And have timepoints. 
    graph = pydot.Dot("my_graph", graph_type='graph', rankdir='TB', label=name)
    gNodeIds = {} # Maintain a dict to with tree nodeId and corresponding Graph node ID
    #usc_nodes = [13, 14, 15]
    colors_list = ['#ffcccc','#ffffb3','#e6ccff','#e6ccff','#e6ccff']
    # First add the nodes in each timepoint
    for tp, nodes in tp_nodes.items():
        if tp == -1:
            gNodeIds[-1] = "Root"
            graph.add_node(pydot.Node("Root", shape="circle", width=1, fontsize=120))
            continue
        colorVal = colors_list[tp]
        #firstNode = True
        # Add nodes to this subgraph for having nodes in same timepoint at a same level
        S = pydot.Subgraph(rank='same')
        for n in nodes:
            nodeId = str(n)+"_"+str(len(id_cells[n]))
            gNodeIds[str(n)] = nodeId

            #if firstNode:
            #    n1 = pydot.Node(nodeId, shape="circle", xlabel="t="+str(tp), fixedsize="true", width=0.7)
                #graph.add_node(pydot.Node(nodeId, shape="circle", xlabel="t="+str(tp), rank="same",fixedsize="true", width=0.7))
           #     S.add_node(n1)
           #     firstNode = False
            if n in usc_nodes: # unobserved subclones are not added in the subgraph
                graph.add_node(pydot.Node(nodeId, shape="circle", style="dotted", width=1, fontsize=120))
            else:
                n2 = pydot.Node(nodeId, shape="circle", color = colorVal, style="filled", width=1, fontsize=120)
                S.add_node(n2)
        graph.add_subgraph(S)

    # Then add the edges and mutations for the nodes
    for node, edge in id_edges.items():
        if node == 0:
            continue
        pId = edge.split('_')[0]
        if pId == '-2':
            continue
        if pId == '0':
            p_gNode = "Root"
        else:
            if pId not in gNodeIds:
                continue
            p_gNode = gNodeIds[pId]

        gNode = gNodeIds[str(node)]
        if id_newMut[node] != set():
            id_newMut[node].sort()
            smutList = []
            # use the following for smaller sample of large SA501
            if sample == "gene":
                lmutList = gene_smallSA501(id_newMut[node])
                smutList = getMutList(lmutList)
                mutations = ";".join(str(i) for i in id_newMut[node])
            if sample == "large":
                smutList = getMutList(id_newMut[node])
                mutations = ";".join(str(i) for i in id_newMut[node])
            if sample == "AML":
                amlMut = AML_mut(id_newMut[node])
                mutations = ";".join(str(i) for i in amlMut)
            else:
                mutations = ";".join(str(i) for i in id_newMut[node])

            #mutGeneList = mutGenes(mutList, id_newMut[node])
            #mutGeneList = set(mutList) & set(id_newMut[node])
            #newMut = set(id_newMut[node]) - set(mutGeneList)
            #print("All mutations ",id_newMut[node])
            #print(mutGeneList," ",newMut)
            #mutations = ";".join(str(i) for i in id_newMut[node])
            #mutations = ";".join(str(i) for i in newMut)
            #if mutGeneList != []:
            #    mutGene = ";".join(str(i) for i in mutGeneList)
            if smutList != []:
                smut = ";".join(str(i) for i in smutList)
            if node in id_backMut:
                backMut = ";".join(str(i) for i in id_backMut[node])
                #graph.add_edge(pydot.Edge(p_gNode, gNode, color="black", label="""<<font color="black">"""+mutations+"""</font><font color="red">"""+backMut+"</font>>"))
                #graph.add_edge(pydot.Edge(p_gNode, gNode, color="black", label="<<font color='black'>"+mutations+"</font><font color='red'>"+backMut+"</font>>"))
                graph.add_edge(pydot.Edge(p_gNode, gNode, color="black", label=mutations, fontsize=120))
                graph.add_edge(pydot.Edge(p_gNode, gNode, color="red", label=backMut, fontcolor="red", fontsize=120))
                if smutList != []:
                    graph.add_edge(pydot.Edge(p_gNode, gNode, color="blue", label=smut, fontcolor="blue", fontsize=120))
            else:
                #if node in id_plMut:
                #    id_plMut[node].sort()
                #    plMut = ";".join(str(i) for i in id_plMut[node])
                #    graph.add_edge(pydot.Edge(p_gNode, gNode, color="blue", label=plMut, fontcolor="blue"))
                #else:
                graph.add_edge(pydot.Edge(p_gNode, gNode, color="black", label=mutations, fontsize=120))
                if smutList != []:
                    graph.add_edge(pydot.Edge(p_gNode, gNode, color="blue", label=smut, fontcolor="blue", fontsize=120))
        else:
            if node in id_backMut:
                backMut = ";".join(str(i) for i in id_backMut[node])
                graph.add_edge(pydot.Edge(p_gNode, gNode, color="black", label=backMut, fontcolor="red", fontsize=120))
            else:
                graph.add_edge(pydot.Edge(p_gNode, gNode, color="black", fontsize=120))

        #if node in id_backMut:
        #    mutations = ";".join(str(i) for i in id_backMut[node])
        #    graph.add_edge(pydot.Edge(p_gNode, gNode, color="red", label="""<<font color="black">abcde</font><font color="red">f</font>>"""))
        
    #graph.add_node(pydot.Node("a", shape="circle", xlabel="X1", fixedsize="true", width=0.3))
    #graph.add_node(pydot.Node("c", shape="circle", fixedsize="true", width=0.3))
    #graph.add_node(pydot.Node("b", shape="circle", xlabel="X2", fixedsize="true", width=0.3))
    #graph.add_node(pydot.Node("d", shape="circle", xlabel="X2", fixedsize="true", width=0.3))
    #graph.add_edge(pydot.Edge("a", "b", color="blue", label="1-12"))
    
    #graph.add_node(pydot.Node("1", label="Root", color="red", style="filled",shape="circle", xlabel="X1", fontsize=10))
    #graph.add_node(pydot.Node("2", label="", color="blue", style="filled",shape="circle", xlabel="X2", fontsize=10))
    #graph.add_edge(pydot.Edge("1", "2", style="dotted", label="1,2,3"))
    #graph.write_png(opFile)
    return graph

#tp_nodes = {-1: [0], 0: [1, 2, 3, 4, 13], 1: [5, 6, 7, 8, 9, 14], 2: [10, 11, 12, 15]}
#id_cells = {0: [], 1: [0, 1, 2, 23], 2: [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], 3: [22], 4: [24, 25, 26], 5: [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40], 6: [41, 42, 44, 45, 46, 47, 48], 7: [43], 8: [49, 50, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62], 9: [51, 52], 10: [63], 11: [64, 66, 67, 68, 69, 70, 71, 72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89], 12: [65, 74, 75], 13: [], 14: [], 15: []}
#id_newMut = {0: set(), 1: {9, 10}, 2: {11, 12, 13, 7}, 3: {0, 17, 2, 5, 6, 15}, 4: {11, 12, 13, 14, 15, 16, 17, 18}, 5: {8, 7}, 6: {0, 1, 3, 4, 5, 7, 11}, 7: set(), 8: set(), 9: {0, 1, 3, 4, 5, 8, 9, 10, 11, 13, 14, 15, 17, 18}, 10: set(), 11: {5, 13, 14, 15, 16, 17, 18, 19}, 12: {0, 1, 3, 11, 13, 15, 16, 17, 18, 19}, 13: {0, 1, 2, 3, 4, 5, 6}, 14: set(), 15: set()}
#id_backMut = {5: [0], 14: [0, 1, 3, 4, 5, 7, 11, 13], 15: [5, 7]}
#id_edges = {0: '-1_0', 1: '13_1', 2: '13_2', 3: '0_3', 4: '13_4', 5: '1_5', 6: '14_6', 7: '14_7', 8: '4_8', 9: '14_9', 10: '15_10', 11: '15_11', 12: '7_12', 13: '0_13', 14: '2_14', 15: '6_15'}

#tp_nodes = {-1: [0], 0: [1, 2, 3, 4, 13], 1: [5, 6, 7, 8, 9], 2: [10, 11, 12]}
#id_cells = {0: [], 1: [0, 1, 2, 23], 2: [3, 4, 5, 6, 7, 8, 9, 10, 11, 12, 13, 14, 15, 16, 17, 18, 19, 20, 21], 3: [22], 4: [24, 25, 26], 5: [27, 28, 29, 30, 31, 32, 33, 34, 35, 36, 37, 38, 39, 40], 6: [41, 42, 44, 45, 46, 47, 48], 7: [43], 8: [49, 50, 53, 54, 55, 56, 57, 58, 59, 60, 61, 62], 9: [51, 52], 10: [63], 11: [64, 66, 67, 68, 69, 70, 71, 72, 73, 76, 77, 78, 79, 80, 81, 82, 83, 84, 85, 86, 87, 88, 89], 12: [65, 74, 75], 13: []}
#id_newMut = {0: [], 1: [9, 10], 2: [11, 12, 13], 3: [], 4: [], 5: [8, 7], 6: [], 7: [], 8: [], 9: [], 10: [], 11: [14, 15, 16, 17, 18, 19], 12: [], 13: [0, 1, 2, 3, 4, 5, 6]}
#id_edges = {0: '-1_0', 1: '13_1', 2: '13_2', 3: '0_3', 4: '13_4', 5: '1_5', 6: '2_6', 7: '2_7', 8: '4_8', 9: '2_9', 10: '6_10', 11: '6_11', 12: '7_12', 13: '0_13'}
#id_backMut = {}

#drawTree(tp_nodes, id_cells, id_newMut, id_backMut, id_edges, "plots/iter"+str(1)+"test.png")
