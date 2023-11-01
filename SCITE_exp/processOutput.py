import argparse
import pandas as pd
from random import randint, choice
import random
import json

''' Get cluster genotype by taking a majority vote for each cell. '''
def vote_cluster_genotypes(cell_gs_list):
    cell_gs_len = len(cell_gs_list)
    cg_len = len(cell_gs_list[0])
    #print(" Inside Voted cluster genotype ")
    #print(cell_gs_len," ",cg_len)
    #print(cell_gs_list[0])
    #print("============================")
    voted_cg = []
    #total0s = 0
    #total1s = 0
    if len(cell_gs_list) == 1:
        for j in range(0,cg_len):
            if cell_gs_list[0][j] == '3':
                k = randint(0, 1)
                voted_cg.append(str(k))
            else:
                voted_cg.append(cell_gs_list[0][j])
    else:    
        for j in range(0,cg_len):
            total0s = 0
            total1s = 0
            for i in range(0,cell_gs_len):
                if cell_gs_list[0][j] == '3':
                    k = randint(0, 1)
                    if k == 0:
                        total0s = total0s+1
                    else:
                        total1s = total1s+1
                if cell_gs_list[i][j] == '0':
                    total0s = total0s+1
                if cell_gs_list[i][j] == '1':
                    total1s = total1s+1
            #print(" Total 0s ",total0s," Total 1s ",total1s)
            if total0s > total1s:
                voted_cg.append('0')
            if total1s >= total0s:
                voted_cg.append('1')
    return voted_cg

def gen_cell_genotype(clones_cells, noOfCells, clonal_genotype, op_gen):
    # Get an array with cells as indices and genotype as values.
    cell_genotype = [clonal_genotype[0]] * int(noOfCells)
    for clones, cells in clones_cells.items():
        clone_genotype = clonal_genotype[int(clones)]
        for cell in cells:
            cell_genotype[int(cell)] = clone_genotype
    #print(clonal_genotype[13])
    #print(" Testing if clones and cells genotypes are same here ")
    #print(clonal_genotype[0] == cell_genotype[315])
    #print(clonal_genotype[0] == cell_genotype[459])
    #print(clonal_genotype[17] == cell_genotype[2])
    #print(clonal_genotype[17] == cell_genotype[28])
    #print(clonal_genotype[8] == cell_genotype[459])
    df = pd.DataFrame(cell_genotype)
    df.to_csv(op_gen, sep='\t')
    return cell_genotype

''' Read the input genotype matrix and generate a genotype matrix with majority voting. '''
def gen_clonal_genotype(input_file, clones_cells):
    f = open(input_file,'r')
    line = f.readline().rstrip('\n')
    D_matrix = []
    while(line != ""):
        genotype = line.split('\t')
        D_matrix.append(genotype)
        line = f.readline().rstrip('\n')
    #print(" D_matrix ")
    #print(D_matrix)
    clonal_genotype = []
    for clone, cells in clones_cells.items():
        cell_gs_list = []
        #print(" Clone ",clone)
        for cell in cells:
            cell_gs_list.append(D_matrix[int(cell)])
        #print(" Cell genotype ",cell_gs_list)
        voted_cg = vote_cluster_genotypes(cell_gs_list)
        clonal_genotype.append(voted_cg)
    #print(len(clonal_genotype))
    return clonal_genotype

''' Create an assignment file to calculate V-measure. '''
def create_assignment_file(clones_cells, noOfCells, op_ass):
    cell_assignment = [-1] * int(noOfCells)
    #print("=============================")
    for clones, cells in clones_cells.items():
        #print("Clone ",clones," cells ",cells)
        for cell in cells:
            #cell_assignment.insert(int(cell), clones)
            cell_assignment[int(cell)] = clones
        #print(cell_assignment)
    #print(cell_assignment)
    f = open(op_ass, "w")
    f.write(' '.join(list(map(str,cell_assignment))))
    f.close()

''' Change the clones to index nos starting from 0 to create an assignment file for calculating V-measure. '''
def change_clone_keys(clones_cells):
    clones_list = list(clones_cells.keys())
    new_clone_keys = {}
    new_clones_cells = {}
    clone_count = 0
    for clone in clones_list:
        new_clone_keys[clone] = clone_count
        clone_count = clone_count+1
    #print(" New clone keys ",new_clone_keys)
    for clone, cells in clones_cells.items():
        clone_no = new_clone_keys.get(clone)
        new_clones_cells[clone_no] = cells
    #print(" New clone cells ",new_clones_cells)
    return new_clones_cells

''' Get the cells belonging to each clone. '''
def map_clones_cells(input_file):
    fileName = open(input_file,'r')
    line = fileName.readline().rstrip('\n')
    clones_cells = {}
    existing_cells = []
    while(line != ""):
        if 'node' not in line and 's' in line:
            clone = (line.split(' -> '))[0]
            cells = ((((line.split(' -> '))[1]).split('s'))[1]).rstrip(';')
            if cells in existing_cells:
                line = fileName.readline().rstrip('\n')
                continue
            if clone in clones_cells:
                cell_list = clones_cells[clone]
                cell_list.append(cells)
                existing_cells.append(cells)
                clones_cells[clone] = cell_list
            else:
                cell_list = []
                cell_list.append(cells)
                existing_cells.append(cells)
                clones_cells[clone] = cell_list
        line = fileName.readline().rstrip('\n')
    #print(" Clones to cells ",clones_cells)
    #print(" Clones ",clones_cells.keys())
    return clones_cells

''' Get the timepoint of the cells. '''
def get_timepoint_cells(cells_file):
    cfile = open(cells_file,'r')
    cline = cfile.readline().rstrip('\n')
    timepoint_cells = {}
    cellCount = 0
    while(cline != ""):
        cline_list = cline.split('\t')
        timepoint = cline_list[0].split('_')[0]
        cells_list = cline_list[1].split(';')
        cellCount = cellCount+len(cells_list)
        prev_timepoint = int(timepoint)-1
        timepoint_key = str(prev_timepoint)+"_"+timepoint

        if timepoint_key in timepoint_cells:
            temp_cells = timepoint_cells[timepoint_key]
            temp_cells.extend(cells_list)
            timepoint_cells[timepoint_key] = temp_cells
        else:
            temp_cells = []
            temp_cells.extend(cells_list)
            timepoint_cells[timepoint_key] = temp_cells

        cline = cfile.readline().rstrip('\n')
    print(" Timepoint cells ",timepoint_cells)
    print(" No. of cells ",cellCount)
    return timepoint_cells, cellCount

''' Get the mutations at each timepoint for evaluation. '''
def get_timepoint_mutations(timepoint_cells, cell_genotype, op_path):
    timepoint_mutations = {}
    for tp, cells in timepoint_cells.items():
        for cell in cells:
            c_genotype = cell_genotype[int(cell)]
            mutations = [i for i in range(len(c_genotype)) if c_genotype[i]=='1']
            #print(" Mutations ",mutations)
            if tp in timepoint_mutations:
                temp_mut = timepoint_mutations[tp]
                temp_mut.extend(mutations)
                timepoint_mutations[tp] = temp_mut
            else:
                temp_mut = []
                temp_mut.extend(mutations)
                timepoint_mutations[tp] = temp_mut
    
    #print(" Timepoint mutations ",timepoint_mutations)
    #print(" Timepoint mutation set ",type(timepoint_mutations))
    w_file = open(op_path+'/tp_mut.json',"w")
    json.dump(timepoint_mutations,w_file)

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input data matrix")
parser.add_argument("-tree", "--tree",dest ="tree", help="Tree output by SCITE")
parser.add_argument("-cells", "--cells",dest ="cells", help="Cells belonging to each timepoint. (Generated in simulation)")
parser.add_argument("-op_ass", "--op_ass",dest ="op_ass", help="Output file to save assignment file")
parser.add_argument("-op_gen", "--op_gen",dest ="op_gen", help="Output file to save consensus genotype")
parser.add_argument("-op_path", "--op_path",dest ="op_path", help="Output path to save the timepoint mutations for evaluation")
args = parser.parse_args()
clone_cells = map_clones_cells(args.tree)
updated_clone_cells = change_clone_keys(clone_cells)
timepoint_cells, noOfCells = get_timepoint_cells(args.cells)

print(" Assignment file name ",args.op_ass)
create_assignment_file(updated_clone_cells, noOfCells, args.op_ass)
clonal_genotype = gen_clonal_genotype(args.input, updated_clone_cells)
cell_genotype = gen_cell_genotype(updated_clone_cells, noOfCells, clonal_genotype, args.op_gen)

# Loop over timepoint_cells and for each cell get the genotype and the mutation based on the index of the genotype array.
get_timepoint_mutations(timepoint_cells, cell_genotype, args.op_path)
# Also, get the consensus genotype based on majority voting.
