import argparse

def inputToDict(graph):
    graph_dict = {}
    file = open(graph,"r")
    line = file.readline().rstrip('\n')
    while(line != ""):
        line_a = line.split('\t')
        nodeId = line_a[0]
        timepoint = line_a[1]
        pIDs = line_a[2]
        muts = line_a[4]
        edge_len = line_a[5]
        graph_dict[nodeId+'_'+timepoint] = pIDs+";"+muts+";"+edge_len
        line = file.readline().rstrip('\n')
    return graph_dict

def get_node_children(graph_dict):
    child_dict = {}
    for key, value in graph_dict.items():
        nodeId = (key.split('_'))[0]
        #print(type((value.split(';'))[0]))
        pIDs = ((value.split(';'))[0]).split(',')
        print(pIDs)
        if not pIDs == []:
            for p in pIDs:
                #if p == '[' or p == ']' or p == ',' or p == ' ':
                #    continue
                print(p)
                if p in child_dict:
                    child_list = child_dict[p]
                    child_list.append(nodeId)
                    child_dict[p] = child_list
                else:
                    child_list = []
                    child_list.append(nodeId)
                    child_dict[p] = child_list
    print(child_dict)

#def trim_unobservedSubclones(graph_dict):


parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Graph with arranged subclones and unobserved subclones.")
args = parser.parse_args()

graph_dict = inputToDict(args.input)
print(graph_dict)
get_node_children(graph_dict)

# Modify the graph dict to update or trim unobserved subclones and edges. 
