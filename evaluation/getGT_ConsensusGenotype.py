import argparse
import pandas as pd

''' Reads the default.SNVcell.csv file. '''
def cells_timepoint(cell_file):
    file = open(cell_file,"r")
    line = file.readline().rstrip('\n')
    tp_node_cells = {}
    timepoint_cells = {}
    while(line != ""):
        line_a = line.split('\t')
        tp_node = line_a[0]
        cells = line_a[1]
        tp_node_cells[tp_node] = cells
        line = file.readline().rstrip('\n')

    for tp_node, cells in tp_node_cells.items():
        timepoint = (tp_node.split('_'))[0]
        cell_list = cells.split(';')
        if timepoint in timepoint_cells:
            c_list = timepoint_cells[timepoint]
            c_list.extend(cell_list)
            timepoint_cells[timepoint] = c_list
        else:
            c_list = []
            c_list.extend(cell_list)
            timepoint_cells[timepoint] = c_list
    #print(timepoint_cells)
    return timepoint_cells

''' Reads the cell_genotype file. '''
def read_cg(cg_file):
    file = open(cg_file,"r")
    line = file.readline().rstrip('\n')
    cellCount = 0
    cell_genotype = {}
    while(line != ""):
        line_a = line.split('\t')
        cell_genotype[cellCount] = line_a
        cellCount = cellCount+1
        line = file.readline().rstrip('\n')
    return cell_genotype

''' Save the cell genotype at each timepoint. '''
def save_cell_genotype(timepoint_cells, cell_genotype, op_file):
    for tp, cells in timepoint_cells.items():
        tp_mutations = []
        for cell in cells:
            cell_mut = cell_genotype[int(cell)]
            #print(cell_mut)
            tp_mutations.append(cell_mut)
        tp_mut_df = pd.DataFrame(tp_mutations)
        print(" Timepoint ",tp," mutations ",tp_mut_df)
        tp_mut_df.to_csv(op_file+"_t"+tp+".csv",sep='\t')

parser = argparse.ArgumentParser()
parser.add_argument("-cg", "--cg",dest ="cg", help="default.G.csv file from simulation")
parser.add_argument("-cell", "--cell",dest ="cell", help="default.SNVcell.csv file from simulation")
parser.add_argument("-op","--op",dest="op", help="Output file to save")
args = parser.parse_args()

timepoint_cells = cells_timepoint(args.cell)
cell_genotype = read_cg(args.cg)
save_cell_genotype(timepoint_cells, cell_genotype, args.op)
