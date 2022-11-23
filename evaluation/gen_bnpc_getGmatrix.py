# This file gets the G matrix with the cells and mutations.
import pandas as pd
import argparse

def get_doublet_cells(cell_clone_input):
    input_list = cell_clone_input.split("/")
    simDataType = input_list[2]
    repNo = input_list[3]
    print("sim_input/"+simDataType+"/"+repNo+"/input_"+simDataType+"_"+repNo+".SNVcell.csv")
    doublet_info_file = "sim_input/"+simDataType+"/"+repNo+"/input_"+simDataType+"_"+repNo+".SNVcell.csv"
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

def doublet_get_cell_cluster(assignment, doublet_cells_list):
    with open(assignment,'r') as f:
        lines = f.readlines()

    assignment_list = ((lines[1].split('\t'))[2]).split(" ")
    cell_cluster = {}
    print(len(assignment_list))
    cell_count = 1
    for i in assignment_list:
        if cell_count in doublet_cells_list:
            continue
        if '\n' in i:
            i = i.replace("\n","")
        cell_cluster[cell_count] = i
        cell_count = cell_count+1
    return cell_cluster,cell_count

''' Cells are the keys and cluster are the values. '''
def get_cell_cluster(assignment):
    with open(assignment,'r') as f:
        lines = f.readlines()

    assignment_list = ((lines[1].split('\t'))[2]).split(" ")
    cell_cluster = {}
    print(len(assignment_list))
    cell_count = 1
    for i in assignment_list:
        if '\n' in i:
            i = i.replace("\n","")
        cell_cluster[cell_count] = i
        cell_count = cell_count+1
    return cell_cluster

''' Get timepoints for each cell from the simulation file. '''
def get_cell_timepoints(celltp_file):
    file = open(celltp_file,"r")
    line = file.readline().rstrip('\n')
    cell_tp = {}
    while(line != ""):
        line_a = line.split('\t')
        #print(line_a)
        tp = (line_a[0].split('_'))[0]
        cells = line_a[1].split(';')
        for cell in cells:
            cell_tp[cell] = tp
        line = file.readline().rstrip('\n')

    tp_Cell = {}
    for cell, tp in cell_tp.items():
        if tp in tp_Cell:
            t_list = tp_Cell[tp]
            t_list.append(cell)
            tp_Cell[tp] = t_list
        else:
            t_list = []
            t_list.append(cell)
            tp_Cell[tp] = t_list
    return tp_Cell

def doublet_get_cell_mutation(mutations, noOfCells, doublet_cells_list):
    cluster_mutation_df = pd.read_csv(mutations, sep='\t', index_col = 0)
    #cluster_mutation_df_T = cluster_mutation_df.T
    print(cluster_mutation_df)
    #input_df = pd.read_csv(input_D, sep='\t', index_col = 0)
    #print(input_df)
    #cells = list(input_df.index) # For doublets remove those cell IDs from the list
    #pos = list(input_df.columns)
    #cluster_mutation_df.index = pos
    cols = [i for i in range(0,noOfCells)]
    cluster_mutation_df.columns = cols
    #mutation_df = cluster_mutation_df.reindex(pos, columns=cells)

    for dCells in doublet_cells_list:
        cluster_mutation_df.drop(int(dCells), inplace=True, axis=1)
    print(cluster_mutation_df)
    return cluster_mutation_df.T

def get_cell_mutation(mutations):
    cluster_mutation_df = pd.read_csv(mutations, sep='\t', index_col = 0)
    print(cluster_mutation_df)
    return cluster_mutation_df.T

def get_timepoints_mutation(tp_cell, mutations, opFile):
    cluster_mutation_df = pd.read_csv(mutations, sep='\t', index_col = 0)
    cluster_mutation_df = cluster_mutation_df.T
    cluster_mutation_list = cluster_mutation_df.values.tolist()
    #print(cluster_mutation_list)
    for tp,cells in tp_cell.items():
        tp_mutations = []
        for cell in cells:
            for i in range(len(cluster_mutation_list)):
                if int(cell) == i:
                    tp_mutations.append(cluster_mutation_list[i])
        tp_mut_df = pd.DataFrame(tp_mutations)
        print(" Timepoint ",tp," mutations ",tp_mut_df)
        tp_mut_df.to_csv(opFile+"_t"+tp+".csv",sep='\t')

parser = argparse.ArgumentParser()
parser.add_argument("-cc", "--cc",dest ="cc", help="Assignment.txt file indicating clusters of cells")
parser.add_argument("-gp", "--gp",dest ="gp", help="Mutations for each cluster")
#parser.add_argument("-input", "--input",dest="input", help="Sample input to bnpc")
parser.add_argument("-celltp", "--celltp",dest="celltp", help="Cell timepoint")
parser.add_argument("-doublet", "--doublet",dest="doublet", help="Doublet data")
parser.add_argument("-op","--op",dest="op", help="Output file to save")
args = parser.parse_args()

if args.doublet == "true":
    doublet_cells = get_doublet_cells(args.cc)
    cell_cluster,cell_count = doublet_get_cell_cluster(args.cc, doublet_cells)
    cell_mutation_df = doublet_get_cell_mutation(args.gp, cell_count-1, doublet_cells)
    cell_mutation_df.to_csv(args.op,sep='\t')
else:
    cell_cluster = get_cell_cluster(args.cc)
    #print(cell_cluster.values())
    cell_mutation_df = get_cell_mutation(args.gp)
    cell_mutation_df.to_csv(args.op+".tsv",sep='\t') # Save the entire consensus genotype.
    timepoints_cell = get_cell_timepoints(args.celltp)
    print(timepoints_cell)
    get_timepoints_mutation(timepoints_cell, args.gp, args.op) # Save the consensus genotype at each timepoint.


