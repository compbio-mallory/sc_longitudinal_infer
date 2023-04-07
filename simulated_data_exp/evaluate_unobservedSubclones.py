import argparse

''' Get the mutations of each node. '''
def get_nodeMut(mut_file):
    file = open(mut_file,"r")
    line = file.readline().rstrip('\n')
    node_mut = {}
    node_mut['0'] = [] # Add the root.
    while(line != ""):
        line_a = line.split('\t')
        node = line_a[0]
        mut = line_a[1].split(';')
        line = file.readline().rstrip('\n')
        node_mut[node] = mut
    return node_mut

''' Get the timepoint of the nodes. '''
def get_tp_node(gtTree):
    file = open(gtTree,"r")
    line = file.readline().rstrip('\n')
    firstLine = True
    node_timepoint = {}
    dead_nodes = []
    while(line != ""):
        if firstLine == True:
            line = file.readline().rstrip('\n')
            firstLine = False
            continue
        line_a = line.split('\t')
        timepoint = line_a[0]
        node = line_a[1]
        leaf = line_a[5]
        node_timepoint[node] = timepoint
        if leaf == '-1':
            dead_nodes.append(node)
        line = file.readline().rstrip('\n')
    return node_timepoint, dead_nodes

def get_unobserved_subclone(node_mut, node_timepoint, dead_nodes):
    usc_mut = {}
    for node, mut in node_mut.items():
        if node in dead_nodes:
            print(" Dead! ",node)
            continue
        node_tp = float(node_timepoint[node])
        if node_tp % 1 == 0:
            continue
        if node_tp in usc_mut:
            t_mut = usc_mut[node_tp]
            t_mut.extend(mut)
            usc_mut[node_tp] = t_mut
        else:
            t_mut = []
            t_mut.extend(mut)
            usc_mut[node_tp] = t_mut # Save the mutations of the unobserved subclones at each timepoint.
    print(" Unobserved subclones ",usc_mut)
    return usc_mut

''' Read the inferred tree and save the mutations, parent node and timepoints. '''
def read_inferred_tree(lgTree):
    file = open(lgTree,'r')
    line = file.readline().rstrip('\n')
    node_mut = {}
    node_timepoint = {}
    child_parent = {}
    while(line != ""):
        line_a = line.split('\t')
        node = line_a[0]
        timepoint = line_a[1]
        node_timepoint[node] = timepoint
        parent = line_a[2]
        child_parent[node] = parent
        if node == '0':
            node_mut[node] = []
        else:
            mut = line_a[4].split(',')
            node_mut[node] = mut
        line = file.readline().rstrip('\n')
    return node_mut, node_timepoint, child_parent

''' Evaluate unobserved subclones. '''
def evaluate_unobservedSubclone(gt_unobserved_subclone, lg_unobserved_subclone):
    for tp, mut in gt_unobserved_subclone.items():
        noOf_gt_mutations = len(mut) # total number of correct mutations for unobserved subclones at a timepoint.
        if tp in lg_unobserved_subclone:
            lg_mut = lg_unobserved_subclone[tp] # mutations from the inferred longitudinal tree.
        else:
            lg_mut = []
        correct_mut = len(set(mut).intersection(set(lg_mut))) # Correct number of mutations inferred by the longitudinal tree algorithm.
        incorrect_mut = noOf_gt_mutations - correct_mut # Incorrect number of mutations inferred by the longitudinal tree algorithm.
        correct_usc = correct_mut / noOf_gt_mutations
        incorrect_usc = incorrect_mut / noOf_gt_mutations
        print(" Timepoint ",tp," Correct usc ",correct_usc," Incorrect usc ",incorrect_usc)

parser = argparse.ArgumentParser()
parser.add_argument("-mut", "--mut",dest ="mut", help="Ground truth mutations")
parser.add_argument("-gtTree", "--gtTree",dest ="gtTree", help="Ground truth tree")
parser.add_argument("-lgTree","--lgTree",dest="lgTree", help="Longitudinal tree inferred from the algorithm")
args = parser.parse_args()

print(" For GT tree ",args.gtTree," longitudinal Tree ",args.lgTree)
gt_node_mut = get_nodeMut(args.mut)
print(" Node mutations ",gt_node_mut)
gt_node_timepoint, dead_nodes = get_tp_node(args.gtTree)
print(" Node timepoint ",gt_node_timepoint)
print(" Dead nodes ",dead_nodes)
gt_unobserved_subclone = get_unobserved_subclone(gt_node_mut, gt_node_timepoint, dead_nodes)

print("============= Longitudinal Tree =================")
lg_node_mut, lg_node_timepoint, lg_child_parent = read_inferred_tree(args.lgTree)
lg_unobserved_subclone = get_unobserved_subclone(lg_node_mut, lg_node_timepoint, [])

evaluate_unobservedSubclone(gt_unobserved_subclone, lg_unobserved_subclone)
