import argparse

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

''' Finds the tree or trees with minimum number of parallel and back mutations. '''
def calculateParallelBackMutations(treeFiles):
    tree_back_parallelMut = {}

    for treef in treeFiles:
        cluster_name = (treef.split('/')[1]).split('_')[0]
        print(" Cluster name ",cluster_name)
        tf = open(treef,'r')
        tf_line = tf.readline().rstrip("\n")
        mut_list = []
        node_mut = {}
        child_parent = {}
        while(tf_line != ""):
            node_id = tf_line.split('\t')[0]
            if node_id == '0': # Skipping the root node here
                tf_line = tf.readline().rstrip("\n")
                continue
            timepoint = tf_line.split('\t')[1]
            pid = tf_line.split('\t')[2]
            mut = (tf_line.split('\t')[4]).split(',')
            mut_list.extend(mut)
            node_mut[node_id] = mut
            child_parent[node_id] = pid
            tf_line = tf.readline().rstrip("\n")
        print(node_mut)
        mut_list = set(mut_list)
        print(" Mutations ",mut_list)

        all_mutations = {}
        parent_child = {}
        new_mutations = {}

        for child, parent in child_parent.items(): # Find the actual mutation set at a timepoint for each node.
            if parent == '0':
                all_mutations[child] = node_mut[child]
                continue
            if "," in parent:
                split_parent = parent.split(",")
                p_mut = node_mut[split_parent[0]]
            else:
                p_mut = node_mut[parent]

            if parent in parent_child:
                child_list = parent_child[parent]
                child_list.append(child)
                parent_child[parent] = child_list
            else:
                parent_child[parent] = [child]
            child_mut = node_mut[child]
            actual_mutation = set(child_mut) - set(p_mut) # Actual mutation at a timepoint
            print(" Actual mutation of ",child," ",list(actual_mutation))
            new_mutations[child] = list(actual_mutation)
            all_mutations[child] = child_mut

        print(" Updated node mutation ",all_mutations)
        parallel_mut = getParallelMutCount(new_mutations)
        back_mut = getBackMutCount(all_mutations, new_mutations, parent_child)
        total_parallel_back = parallel_mut + back_mut
        tree_back_parallelMut[treef] = total_parallel_back
    
    treeFile_ForMin_mutations = min(tree_back_parallelMut, key=tree_back_parallelMut.get)
    print(" Best Tree file ",treeFile_ForMin_mutations)

parser = argparse.ArgumentParser()
parser.add_argument("-tree", "--tree",dest ="tree", nargs="*", help="Final Longitudinal trees")
args = parser.parse_args()

print(" List of tree files ",args.tree)
calculateParallelBackMutations(args.tree)
