import argparse
import sys
from sklearn.metrics.cluster import v_measure_score
from usage import bnpc_assignment, scg_assignment, scclone_assignment, robustclone_assignment, scite_assignment, siclonefit_assignment

# make clu_supp each row composed of the cells in this cluster separated by ;
def read_clu_results(file):
    ret = []
    total = 0
    f = open(file, "r")
    line = f.readline().rstrip('\n')
    while(line != ""):
        ret.append(line)
        cell_num = len(line.split(";"))
        total = total + cell_num
        line = f.readline().rstrip('\n')
    f.close()
    return ret, total

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
# if with_header is false, then no extra row or column
def gen_GT_clu(gt_file, with_header):
    ret = []
    file = open(gt_file, "r")
    line = file.readline().rstrip('\n')
    if line != "" and with_header:
        # data starts from the second row
        line = file.readline().rstrip('\n')
    elif line == "":
        return ret
    l = 0
    gt_dict = {}
    while(line != ""):
        line_a = line.split('\t')
        gt = ""
        if with_header:
            gt = "".join(line_a[1:])
        else:
            gt = "".join(line_a)

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

    return ret

def get_doublet_cells(doublet_info_file):
    with open(doublet_info_file,'r') as fp:
        lines = fp.readlines()
        #print(lines)

    doublet_cells_list = []
    for line in lines:
        first_col = line.split('\t')
        if ';' in first_col[0]:
            cells = first_col[1]
            if ';' in cells:
                cells_list = (cells.strip()).split(';')
                doublet_cells_list.extend(cells_list)
            else:
                doublet_cells_list.append(cells.strip())
    print(" Doublet cells ",doublet_cells_list)
    return doublet_cells_list

def remove_doublet_CG_GT(clu_IDs, gt_IDs, doublet_cells_list):
    print(len(clu_IDs))
    print(len(gt_IDs))
    print(doublet_cells_list)
    count = 0
    for d_cells in doublet_cells_list:
        #del clu_IDs[int(d_cells)]
        #del gt_IDs[int(d_cells)]
        dcell_index = int(d_cells) - count
        clu_IDs.pop(dcell_index)
        gt_IDs.pop(dcell_index)
        count = count+1
    return clu_IDs, gt_IDs

if len(sys.argv) <= 1:
    print("""
    Evaluate clustering results. 
    Usage: python evaluate.py
        -i (--input-file)   Provides the inferred clustering results. (NA) 
        -G (--gt-file)      Provides the ground truth genotype file. (NA)
        -H (--with-header)  Whether the ground truth genotype file has an extra row and column as the header. (False)
        -v (--v-measure)    Select v-measure. (True)
        -d (--doublet)      Indicates presence or absence of doublets. (False)
        -df (--doublet-file)The file mentioning the cells having doublets. (NA)
        """)

    sys.exit(0)

# TODO: turn input file into a string such as bnpc:file_location;scg:file_location, and return the v-measure for only these files. 

parser = argparse.ArgumentParser(description='Evaluate clustering results. ')
parser.add_argument('-i', '--input-files', default="NA")
parser.add_argument('-G', '--gt-file', default="NA")
parser.add_argument('-d', '--doublet', default="NA")
parser.add_argument('-df', '--doubletFile', default="NA")
parser.add_argument('-H', '--with-header', action='store_true')
parser.add_argument('-v', '--v-measure', action='store_true')

args = parser.parse_args()
inputs = args.input_files.split(";")
print(" ================================= ")
print(inputs)
print(" ================================= ")
gt_f = args.gt_file
# whether the ground truth file has one more row and column as the header or not (always assume it is cell by mutation no matter what)
with_header = args.with_header
doublet = args.doublet
print(" Doublet ",doublet)
#print(str(with_header))
if_v_measure = args.v_measure


# ground truth
gt_IDs = gen_GT_clu(gt_f, with_header)

# dissect inputs
for i in inputs:
    method, result_file = i.split(":")
    clu_IDs = []
    num_cells = 0

    if method == "sccluster":
        # inferred clustering results
        clu_supp, num_cells = read_clu_results(result_file)
        clu_IDs = turn_to_clusterIDs(clu_supp, num_cells) 

    elif method == "scg":
        clu_IDs = scg_assignment(result_file)
        num_cells = len(clu_IDs) 

    elif method == "scclone":
        clu_IDs = scclone_assignment(result_file, with_header)
        num_cells = len(clu_IDs) 

    elif method == "bnpc":
        clu_IDs = bnpc_assignment(result_file)
        num_cells = len(clu_IDs)

    elif method == "robustclone":
        clu_IDs = robustclone_assignment(result_file)
        num_cells = len(clu_IDs)

    elif method == "scite":
        clu_IDs = scite_assignment(result_file)
        num_cells = len(clu_IDs)

    elif method == "siclonefit":
        clu_IDs = siclonefit_assignment(result_file)
        num_cells = len(clu_IDs)

    if if_v_measure:
        if doublet == "true":
            doublet_cells_list = get_doublet_cells(args.doubletFile)
            clu_IDs, gt_IDs = remove_doublet_CG_GT(clu_IDs, gt_IDs, doublet_cells_list)
        print("Method is " + method + ". Number of cells from the inferred result: " + str(len(clu_IDs)) + ", versus that of ground truth: " + str(len(gt_IDs)))
        print(" Cluster IDs ",clu_IDs)
        print(" GT cluster IDs ",gt_IDs)
        v_measure = V_measure(clu_IDs, gt_IDs)
        print("V measure results for " + method + ": " + str(v_measure))
