import pydot
import argparse
import random
import numpy as np

#graph = pydot.Dot("my_graph", graph_type="graph")

# Include same color to nodes in the same timepoint
# Add nodes at same timpoint same colors
#my_node = pydot.Node("a", label="Root", xlabel="t=1") # Instead of foo use the node numbers in the graph
#graph.add_node(my_node)
# Or, without using an intermediate variable:
#graph.add_node(pydot.Node("b", shape="circle", xlabel="t=2", fixedsize="true", width=0.3)) # use the xlabels to set the timeline. Change node sizes according to # of cells.

# Add edges
#my_edge = pydot.Edge("a", "b", color="blue", label="1-12") # Here label can be used as mutations
#graph.add_edge(my_edge)

#graph.add_node(pydot.Node("c", shape="circle", xlabel="t=3"))
#graph.add_node(pydot.Node("d", shape="circle"))

# Or, without using an intermediate variable:
#graph.add_edge(pydot.Edge("b", "c", color="blue"))

#graph.add_edge(pydot.Edge("b", "d", style="dotted"))

#graph.write_png("output.png")

#def checkConsecutive(l):
#    return sorted(l) == list(range(min(l), max(l)+1))

# Get the graph to visualize the tree and save the image
def getGraph(node_mut_count, edge_dict, persistent_nodes, timepoint_node, node_cells, opFile):
    graph = pydot.Dot("my_graph", graph_type="graph")
    colors_list = ['#ffcccc','#ffffb3','#e6ccff','#ccf2ff','#ccffcc','#ccffcc']
    check_color_index = [] # maintain a list to check same color is not assigned to every node

    # First just add all the nodes and their respective timepoints
    for tp, nodes in timepoint_node.items():
        nc = 1 # Just a variable to keep count of the nodes added
        time_color = random.randint(1, 5)
        while time_color in check_color_index: # Continue checking until a unique color is selected
            time_color = random.randint(1, 5) # Generate random color for each timepoint
        check_color_index.append(time_color)
        print(" Time color ",colors_list[time_color])

        for node in nodes:
            if tp == 1.0: # Root node will always be one so this will be different.
                graph.add_node(pydot.Node(node, label="Root", color=colors_list[0], style="filled", shape="circle", xlabel="t="+str(tp), fontsize=10))
            else:
                if nc == 1:
                    if tp.is_integer() and node != '0': # Only vary the cells in each integer timepoint as unobserved subclones doesn't have cells sequenced
                        cells_node = node_cells[int(node)] # cell distribution % of the node
                        print(" Node ",node," cells ",cells_node)
                        graph.add_node(pydot.Node(node, label=round(cells_node,2), color=colors_list[time_color], style="filled", shape="circle", xlabel="t="+str(tp), fixedsize="true", width=round(cells_node,2), fontsize=10))
                    else:
                        graph.add_node(pydot.Node(node, label="0", labelloc="b", color=colors_list[time_color], style="filled", shape="circle", xlabel="t="+str(tp), fontsize=10))
                    nc=nc+1
                else:
                    if tp.is_integer():
                        cells_node = node_cells[int(node)] # cell distribution % of the node
                        print(" Node ",node," cells ",cells_node)
                        graph.add_node(pydot.Node(node, label=round(cells_node,2), color=colors_list[time_color], style="filled", shape="circle", fixedsize="true", width=round(cells_node,2), fontsize=10))
                    else:
                        graph.add_node(pydot.Node(node, label="0", labelloc="b", color=colors_list[time_color], style="filled", shape="circle", fontsize=10))

    # Then add the edges and other things
    for tp, nodes in timepoint_node.items():
        for node in nodes:
            if node not in edge_dict:
                continue
            child_nodes = edge_dict[node]
            for child in child_nodes:
                if child in persistent_nodes:
                    graph.add_edge(pydot.Edge(node, child, style="dotted"))
                else:
                    mut_count = node_mut_count[child]
                    graph.add_edge(pydot.Edge(node, child, color="blue", headlabel=mut_count, labeldistance="2", labelfontcolor="red", fontsize=10))

    graph.write_png(opFile)

# Get the difference in node mutations first and then the range.
def getMutationCount(node_mut, edge_dict):
    node_mut_count = {}
    for pID, child_list in edge_dict.items():
        if pID == '-1':
            continue
        pID_mut = set(node_mut[pID])
        for child in child_list:
            if child not in node_mut:
                continue
            child_mut = set(node_mut[child])
            diff_mut = child_mut - pID_mut
            node_mut[child] = list(diff_mut)
    print(" Updated node mutation ",node_mut)
    # Get the count of new mutations. It is hard to visualize the mutations if there are more no. and not consecutive. 

    for node, mut in node_mut.items():
        node_mut_count[node] = len(mut)
    return node_mut_count

# Sort the nodes according to the timepoints to assign the nodes in an order
def sortTimepoint(timepoint_node):
    timepoint = list(timepoint_node.keys())
    timepoint.sort()
    sorted_dict = {i: timepoint_node[i] for i in timepoint}
    return sorted_dict

# Get the distribution of cells for each node at each timepoint
def getNoOfCellsForNodes(cloneNode_dict, cloneCells_dict): 
    node_cells = {} # Save the % of cells at each timepoint
    timepoint_cells = {} # Save the no. of cells at each timepoint
    for clone, cells in cloneCells_dict.items():
        tp = clone.split('_')[1]
        if tp in timepoint_cells:
            cellCount = timepoint_cells[tp]
            cellCount = cellCount+cells
            timepoint_cells[tp] = cellCount
        else:
            timepoint_cells[tp] = cells

    print(" Timepoint cells ",timepoint_cells)
    for clone, node in cloneNode_dict.items():
        tp = clone.split('_')[1]
        tp_total_cells = timepoint_cells[tp]
        clone_cells = cloneCells_dict[clone]
        #print(" Total cells ",tp_total_cells," clone_cells ",clone_cells)
        node_cells[node] = float(clone_cells/tp_total_cells)
    print(" Distribution of cells in nodes ",node_cells)
    return node_cells

def readTreeFile(tree_file):
    # From the file get the list of nodes
    # Create a dictionary of parent(key) -> child(value) to get the list of edges
    # Choose first node and create a dictionary to get the timpoint for that node. Node will be the key.
    # Create a list for persistent nodes or nodes with no mutations. Need to include the style based on this.
    # Create a dictionary for the mutations with node as key and mutations included as a range.
    treef = open(tree_file,'r')
    treef_line = treef.readline().rstrip('\n')

    edge_dict = {} # parent is the key and list values are its children.
    timepoint_node = {} # This dict is important to arrange the nodes according to timepoints.
    persistent_nodes = [] # Includes nodes with no new mutations.
    node_mut = {} # Contains mutations of the nodes.
    node_mut_count = {} # Contains the new mutation count for each node rather the edge length.

    while(treef_line != ""):
        line = treef_line.split("\t")
        nodeId = line[0]
        timepoint = line[1]
        mut_count = line[5]

        # Add the nodes to the timepoint and later sort the timepoints to add the nodes.
        if float(timepoint) in timepoint_node:
            n_list = timepoint_node[float(timepoint)]
            n_list.append(nodeId)
            timepoint_node[float(timepoint)] = n_list
        else:
            n_list = []
            n_list.append(nodeId)
            timepoint_node[float(timepoint)] = n_list

        pId = line[2]
        if pId in edge_dict: # Add the child nodes of the parent to get the edges.
            c_list = edge_dict[pId]
            c_list.append(nodeId)
            edge_dict[pId] = c_list
        else:
            c_list = []
            c_list.append(nodeId)
            edge_dict[pId] = c_list

        mut_count = line[5]
        if mut_count == '0':
            persistent_nodes.append(nodeId)
        if mut_count != '':
            node_mut_count[nodeId] = int(mut_count)

        if nodeId not in persistent_nodes: # Not interested to add mutations of nodes with no new mutations.
            mutations = line[4].split(',')
            node_mut[nodeId] = mutations

        treef_line = treef.readline().rstrip('\n')

    timepoint_node = sortTimepoint(timepoint_node) # Return the sorted timepoint for easy assignment of the nodes in the graph
    print(" Edge dict ",edge_dict)
    print(" Persistent nodes ",persistent_nodes)
    print(" Node mutations ",node_mut)
    print(" Timepoint node ",timepoint_node)
    return edge_dict, persistent_nodes, node_mut, node_mut_count, timepoint_node

parser = argparse.ArgumentParser()
parser.add_argument("-tree", "--tree",dest="tree", help="Longitudinal Tree file")
parser.add_argument("-pngFile", "--pngFile",dest="pngFile", help="The png file name to save.")
parser.add_argument("-cloneNode", "--cloneNode",dest="cloneNode", help="Clone to nodes dictionary.")
parser.add_argument("-cloneCells", "--cloneCells",dest="cloneCells", help="Clone to cells dictionary.")
args = parser.parse_args()

edge_dict, persistent_nodes, node_mut, node_mut_count, timepoint_node = readTreeFile(args.tree)
print("------------------------")
#node_mut_count = getMutationCount(node_mut, edge_dict)
print(" Node mutation count ",node_mut_count)
cloneNode_dict = np.load(args.cloneNode,allow_pickle='TRUE').item()
cloneCells_dict = np.load(args.cloneCells,allow_pickle='TRUE').item()
print(" Clone node ",cloneNode_dict)
print(" Clone dict ",cloneCells_dict)
node_cells = getNoOfCellsForNodes(cloneNode_dict, cloneCells_dict)
getGraph(node_mut_count, edge_dict, persistent_nodes, timepoint_node, node_cells, args.pngFile)
