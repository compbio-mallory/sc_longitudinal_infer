import argparse
import numpy as np
import json

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
    child_parent = {}
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
        if leaf == '-1' and float(timepoint) % 1 != 0: # Dead nodes where cells are not sequenced.
            dead_nodes.append(node)
        pID = line_a[2]
        node_timepoint[node] = timepoint
        child_parent[node] = pID
        line = file.readline().rstrip('\n')
    return node_timepoint, child_parent, dead_nodes

''' Update node mutations where if the mutations are same as its parent then its edge will have 0 mutations. '''
def update_node_mut(gt_node_mut,child_parent):
    for child, parent in child_parent.items():
        if int(child) == 0:
            continue
        child_mut = gt_node_mut[child]
        parent_mut = gt_node_mut[parent]
        if len(child_mut) == 1 and child_mut[0] == '': # Check for empty mutations which shouldn't be the case here.
            child_mut = []
            child_mut.extend(parent_mut)
        else:
            child_mut.extend(parent_mut)
        gt_node_mut[child] = child_mut
        #if set(child_mut) == set(parent_mut):
            #print("Child ",child," Parent ",parent)
            #gt_node_mut[child] = []
        #if len(child_mut) > len(parent_mut):
        #    updated_mut = set(child_mut) - set(parent_mut)
        #    gt_node_mut[child] = list(updated_mut)
    return gt_node_mut

''' Map the timepoints whose mutations should be merged. '''
def map_timepoints(timepoints):
    max_timepoint = max(timepoints) # Ignore adding the max timepoint.
    print(" Max timepoint ",max_timepoint)
    map_tp = {}
    for tp in timepoints:
        if tp in map_tp or tp == max_timepoint:
            continue
        if tp % 1 == 0: # checks only for integer timepoints
            mid_tp = (tp + tp + 1)/2
            map_tp[tp] = [tp,mid_tp,tp+1]
    #print(" Map timepoints ",map_tp)
    return map_tp

''' Save mutations for each timepoint. Key of dictionary tp_mutations is timepoint and key = 1 indicates mutations of t=1 and 1.5 till t=2. '''
def get_tp_mutations(node_mut, node_timepoint, dead_nodes):
    tp_mutations = {}
    for node, mut in node_mut.items():
        if node in dead_nodes:
            continue
        node_tp = float(node_timepoint[node])
        #print("Node ",node," Mut ",mut)
        #if not node_tp%1 == 0:
        #    print(node_tp)
        #    node_tp = node_tp - 0.5
        key = node_tp
        if key in tp_mutations:
            mut_list = tp_mutations[key]
            mut_list.extend(mut)
            tp_mutations[key] = mut_list
        else:
            mut_list = []
            mut_list.extend(mut)
            tp_mutations[key] = mut_list
    #print(" Time point mutations ",tp_mutations.keys())
    map_tp = map_timepoints(set(tp_mutations.keys()))

    union_tp_mutations = {}
    for mt, tp in map_tp.items():
        max_tp = max(tp)
        dict_key = str(int(mt))+"_"+str(int(max_tp))
        for t1 in tp:
            if mt == t1 or t1 not in tp_mutations: # same timepoints mutations are not added to the union dict.
                continue
            tp_mut = tp_mutations[t1] 
            if dict_key in union_tp_mutations:
                temp_list = union_tp_mutations[dict_key]
                temp_list.extend(tp_mut)
                union_tp_mutations[dict_key] = temp_list
            else:
                temp_list = []
                temp_list.extend(tp_mut)
                union_tp_mutations[dict_key] = temp_list

    #print(" Union tp mutations ",union_tp_mutations)
    return union_tp_mutations

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

''' Update mutation list to keep only the new mutations at each timepoint. '''
#def get_updated_mutations_forCG(timepoint_mutation):
#    for tp, mut in timepoint_mutation.items():
#        timepoint = int((tp.split("_"))[1])
#        prev_timepoint = timepoint-1
#        prev_prev_timepoint = timepoint-2
#        prev_key = str(prev_prev_timepoint)+"_"+str(prev_timepoint)
#        if prev_key in timepoint_mutation:
#            prev_time_mut = timepoint_mutation[prev_key]
#        else:
#            continue
        #print(" Previous time mut ",mut," Timepoint mut ",prev_time_mut)
#        if set(mut).issubset(set(prev_time_mut)):
#            timepoint_mutation[tp] = []
#    return timepoint_mutation

''' Get the mutations at each timepoint from BnpC's consensus genotype. '''
def get_mutations_fromCG(cgFile, clones_nodes):
    file = open(cgFile,'r')
    line = file.readline().rstrip('\n')
    timepoint_mutation = {}
    clones_list = list(clones_nodes.keys())
    #clones_list = {'3_2': 1, '4_2': 2, '2_2': 3, '0_2': 4}
    print(" Clones ",clones_list)
    while(line != ""):
        line_a = line.split("\t")
        clone = line_a[0]
        #print(" BnpC clone ",clone)
        if clone not in clones_list:
            line = file.readline().rstrip('\n')
            continue
        timepoint = int((line_a[0].split("_"))[1])
        genotype = line_a[1:]
        #print("Genotype ",genotype)
        key = str(timepoint-1)+"_"+str(timepoint)
        #print(" Timepoint ",timepoint," Genotype ",genotype)
        mut_arr = []
        for j in range(len(genotype)):
            if int(genotype[j]) == 1:
                mut_arr.append(j)

        if key in timepoint_mutation:
            t_mut = timepoint_mutation[key]
            t_mut.extend(mut_arr)
            timepoint_mutation[key] = t_mut
        else:
            timepoint_mutation[key] = mut_arr
        line = file.readline().rstrip('\n')
    #print(" Timepoint_mutation ",timepoint_mutation)
    #timepoint_mutation = get_updated_mutations_forCG(timepoint_mutation)
    return timepoint_mutation

''' Calculate the precision and recall of the whole tree. '''
def calculate_precision_recall(gt_tp_mutations, lg_tp_mutations):
    TP = 0
    FP = 0
    FN = 0
    for timepoint, mutations in gt_tp_mutations.items():
        print(" Timepoint ",timepoint)
        mutations = [int(i) for i in mutations]
        set_mutations = set(mutations)
        print(" GT mutations ",set_mutations)
        #print(" Int mutations ",mutations)
        if timepoint in lg_tp_mutations:
            lg_mut = lg_tp_mutations[timepoint]
        lg_mut = [int(i) for i in lg_mut]
        set_lg_mut = set(lg_mut)
        print(" LG mutations ",set_lg_mut)
        TP = TP+len(set_mutations.intersection(set_lg_mut))
        FP = FP+len(set_lg_mut - set_mutations)
        FN = FN+len(set_mutations - set_lg_mut)
        print(" TP ",TP," FP ",FP," FN ",FN)
        
    if TP + FP == 0:
        precision = 0
    else:
        precision = TP / (TP + FP)
    if TP + FN == 0:
        recall = 0
    else:
        recall = TP / (TP + FN)
    print(" Total Precision ",precision)
    print(" Total Recall ",recall)

''' Compare inferred tree mutations with ground truth mutations and calculate precision, recall. '''
def tp_calculate_precision_recall(gt_tp_mutations, lg_tp_mutations):
    timepoint_precision_recall = {}
    for timepoint, mutations in gt_tp_mutations.items():
        print(" Timepoint ",timepoint)
        mutations = [int(i) for i in mutations]
        set_mutations = set(mutations)
        print(" GT mutations ",set_mutations)
        #print(" LG TP mutations ",lg_tp_mutations)
        if timepoint in lg_tp_mutations:
            lg_mut = lg_tp_mutations[timepoint]
        print(" LG mut ",set(lg_mut))
        lg_mut = [int(i) for i in lg_mut]
        set_lg_mut = set(lg_mut)
        print(" LG mutations ",set_lg_mut)
        TP = len(set_mutations.intersection(set_lg_mut))
        FP = len(set_lg_mut - set_mutations)
        FN = len(set_mutations - set_lg_mut)
        print(" TP ",TP," FP ",FP," FN ",FN)
        if TP + FP == 0:
            precision = 0
        else:
            precision = TP / (TP + FP)
        if TP + FN == 0:
            recall = 0
        else:
            recall = TP / (TP + FN)
        print(" Timepoint ",timepoint," Precision ",precision," Recall ",recall)
        timepoint_precision_recall[timepoint] = str(precision)+","+str(recall)
    return timepoint_precision_recall

''' Read the mutations at each timpoint given by LACE '''
def getLACE_mutations(lace_tp_mutFile):
    f = open(lace_tp_mutFile,'r')
    lace_tp_mut = json.load(f)
    print(" LACE mutations ",lace_tp_mut)
    return lace_tp_mut

parser = argparse.ArgumentParser()
parser.add_argument("-mut", "--mut",dest ="mut", help="Ground truth mutations")
parser.add_argument("-gtTree", "--gtTree",dest ="gtTree", help="Ground truth tree")
parser.add_argument("-lgTree","--lgTree",dest="lgTree", help="Longitudinal tree inferred from the algorithm")
parser.add_argument("-cg","--cg",dest="cg", help="BnpC's inferred Gprime or consensus genotype")
parser.add_argument("-clone","--clone",dest="clone", help="Clones mapped to the inferred longitudinal tree")
parser.add_argument("-lace","--lace",dest="lace", help="LACE file is passed or not for evaluation. Value is true or false.")
parser.add_argument("-laceFile","--laceFile",dest="laceFile", help="LACE file having the mutations in each timepoint")
args = parser.parse_args()

gt_node_mut = get_nodeMut(args.mut)
#print(" Node mutations ",gt_node_mut)
gt_node_timepoint, gt_child_parent, dead_nodes = get_tp_node(args.gtTree)
#print(" Node timepoint ",gt_node_timepoint)
#print(" Child parent ",gt_child_parent)
updated_node_mut = update_node_mut(gt_node_mut, gt_child_parent)
gt_tp_mutations = get_tp_mutations(updated_node_mut, gt_node_timepoint, dead_nodes) # use this if you want edge length to be 0.
#gt_tp_mutations = get_tp_mutations(gt_node_mut, gt_node_timepoint)
print(" GT Timepoint mutations ",gt_tp_mutations)

lg_node_mut, lg_node_timepoint, lg_child_parent = read_inferred_tree(args.lgTree)
print("======== Longitudinal Tree ===========")
print(" Node mut ",lg_node_mut)
#print(" Node timepoint ",lg_node_timepoint)
#print(" Node child parent ",lg_child_parent)
updated_lg_node_mut = update_node_mut(lg_node_mut, lg_child_parent)
print(" Updated lg mut ",updated_lg_node_mut)
lg_tp_mutations = get_tp_mutations(updated_lg_node_mut, lg_node_timepoint, [])
print(" LG mutations ",lg_tp_mutations)

#clones_nodes = np.load(args.clone,allow_pickle='TRUE').item() # Clones mapped to the inferred tree nodes.
#print(" Clone to nodes mapping ",clones_nodes)
#bnpc_timepoint_mutation = get_mutations_fromCG(args.cg, clones_nodes) # Evaluating against BnpC's consensus genotype. Dropping clones that have no supporting cells.
#print(" BnpC Timepoint mutations ",bnpc_timepoint_mutation)
#timepoint_precision_recall = calculate_precision_recall(gt_tp_mutations, bnpc_timepoint_mutation)

print(" For inferred tree ",args.lgTree)
timepoint_precision_recall = tp_calculate_precision_recall(gt_tp_mutations,lg_tp_mutations)
calculate_precision_recall(gt_tp_mutations,lg_tp_mutations)
print("=============================")
if args.lace == "true":
    lace_tp_mut = getLACE_mutations(args.laceFile)
    lace_timepoint_precision_recall = tp_calculate_precision_recall(gt_tp_mutations,lace_tp_mut)
    calculate_precision_recall(gt_tp_mutations,lace_tp_mut)
