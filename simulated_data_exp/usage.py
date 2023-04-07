import copy
from sklearn.metrics.cluster import v_measure_score
import numpy

# initialize a matrix that has r rows and c columns
def init_m(r, c, value):
    m = []
    for i in range(r):
        m.append(copy.deepcopy([value]*c))
    return m

def avg(a):
    s = 0
    for i in a:
        s = s + i
    return s / float(len(a))

# turn the cluster results (an array each entry represents a cluster and has the cells belonging to it separated by ;) into an array in which each entry represents the ID of the cell's cluster. n is the number of cells in total. 
def turn_to_clusterIDs(a, n):
    ret = [0]*n
    for k in range(len(a)):
        cells = a[k].split(";")
        for c in cells:
            ret[int(c)] = k
    return ret

# a and b are two arrays, each entry denoting the cluster ID the cell belongs to. The order of the cluster does not matter, but the homogeneity and compleness matter to V measure on how close a is to b. 
def V_measure(a, b):
    return v_measure_score(a, b)

# given a ground truth file with the first row and first column listing the mutation and cell names, the other entries showing the ground truth genotype of the cell on the mutation, output the clustering, i.e., the label of each cell in an array.
# also return the ground truth genotype profile in the order of the cells (a matrix cell by mutation)
def gen_GT_clu(gt_file):
    # the cluster the cells belong to
    ret = []
    # ground truth genotype matrix
    gt_arr = []
    l = 0
    # from genotype profile to the cells belong to this cluster
    gt_dict = {}
    if gt_file[-4:] == ".npy":
        gt_arr = numpy.load(gt_file)
        for i in range(len(gt_arr)):
            gt = "".join(str(gt_arr[i]))
            if gt not in gt_dict.keys():
                gt_dict[gt] = str(i)
            else:
                gt_dict[gt] = gt_dict[gt] + ";" + str(i)
        l = len(gt_arr)
    else:
        file = open(gt_file, "r")
        line = file.readline().rstrip('\n')
        if line != "":
            # data starts from the second row
            line = file.readline().rstrip('\n')
        else:
            return ret
        while(line != ""):
            line_a = line.split('\t')
            gt_arr.append([int(x) for x in line_a[1:]])
            gt = "".join(line_a[1:])
            if gt not in gt_dict.keys():
                gt_dict[gt] = str(l)
            else:
                gt_dict[gt] = gt_dict[gt] + ";" + str(l)
            l = l + 1
            
            line = file.readline().rstrip('\n')
        file.close()

    ret = [0]*l

    # turn the dict to the ret array
    # k is the cluster genotype
    # cluID is the cluster ID
    cluID = 0
    for k in gt_dict.keys():
        cellIDs = gt_dict[k].split(";")
        for i in cellIDs:
            ret[int(i)] = cluID
        cluID = cluID + 1

    return ret, gt_dict, gt_arr

# given a bnpc assignment file that has the assignment of the cells on row 1 column 2 (0-based) separated by space, output the array that has the assignment. 
def bnpc_assignment(ass_file):
    ret = []
    file = open(ass_file, "r")
    line = file.readline().rstrip('\n')
    if line != "":
        # data starts from the second row
        line = file.readline().rstrip('\n')
    else:
        return ret

    file.close()

    # get the assignment
    line_a = line.split('\t')
    ass = line.split('\t')[2]
    ret = ass.split(' ')
    return ret

# given a scclone assignment file that has assignment of the cells separated by space on row 0 column 0, output the array that has the assignment.
def scclone_assignment(ass_file, wheader):
    ret = []
    file = open(ass_file, "r")
    line = file.readline().rstrip('\n')
    ret = line.split(' ')
    if wheader:
        ret = ret[1:]
    return ret

# given a scite assignment file that has assignment of cells separated by space.
def scite_assignment(ass_file):
    ret = []
    file = open(ass_file, "r")
    line = file.readline().rstrip('\n')
    ret = line.split(' ')
    print(" Assignment file ",ret)
    return ret

# given a robustclone assignment file that has assignment of the cells to each clone where each row represents a clone. Output an array where indices are cells and value is the clone.
def robustclone_assignment(ass_file):
    ret = []
    file = open(ass_file, "r")
    line = file.readline().rstrip('\n')
    ret = line.split(' ')
    print(" Assignment file ",ret)
    return ret

# given a siclonefit assignment file that has assignment of cells separated by space.
def siclonefit_assignment(ass_file):
    ret = []
    file = open(ass_file, "r")
    line = file.readline().rstrip('\n')
    ret = line.split(' ')
    print(" Assignment file ",ret)
    return ret

# given a posterior assignment file whose rows are the cells and columns are the clusters, each entry (starting from row 1 and column 1, 0-based) has the posterior probability of the cell belonging to the cluster, return an array with the assignment of each cell in the same order of the rows. 
def scg_assignment(ass_file):
    ret = []
    file = open(ass_file, "r")
    line = file.readline().rstrip('\n')
    if line != "":
        # data starts from the second row
        line = file.readline().rstrip('\n')
    else:
        return ret

    while(line != ""):
        line_a = line.split('\t')

        post = line_a[1:]
        max_p = 0
        cluster = -1
        # select the cluster that has the highest posterior probability
        for i in range(len(post)):
            if float(post[i]) > max_p:
                max_p = float(post[i])
                cluster = i
        ret.append(cluster)

        line = file.readline().rstrip('\n')

    file.close()

    return ret
