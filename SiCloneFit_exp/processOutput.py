import pandas as pd
import argparse
import json

''' Read SiCloneFit's output genotype. '''
def process_input(input_file, output_file):
    file = open(input_file,"r")
    line = file.readline().rstrip('\n')
    modified_lines = []
    while(line != ""):
        line = line.split(" ")
        new_lines = line[1:]
        new_lines_1 = "\t".join(new_lines)
        modified_lines.append(new_lines_1)
        line = file.readline().rstrip('\n')
    with open(output_file,"w+") as f:
        for l in modified_lines:
            f.write(l)
            f.write("\n")

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
    #print(" Timepoint cells ",timepoint_cells)
    #print(" No. of cells ",cellCount)
    return timepoint_cells, cellCount

''' Get the mutations at each timepoint for evaluation. '''
def get_timepoint_mutations(timepoint_cells, cell_genotype, op_path):
    timepoint_mutations = {}
    for tp, cells in timepoint_cells.items():
        for cell in cells:
            c_genotype = cell_genotype[int(cell)]
            mutations = [i for i in range(len(c_genotype)) if c_genotype[i]==1]
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
parser.add_argument("-cells", "--cells",dest ="cells", help="Cells belonging to each timepoint. (Generated in simulation)")
parser.add_argument("-output", "--output",dest ="output", help="Output data matrix")
parser.add_argument("-op_path", "--op_path",dest ="op_path", help="Output path to save the timepoint mutations for evaluation")
args = parser.parse_args()

process_input(args.input, args.output)
input_df = pd.read_csv(args.output, sep='\t', header=None)
#print(input_df)
transposed_df = input_df.T
transposed_df.to_csv(args.output, sep='\t', index=False, header=False)
# Need to convert the transposed_df into an array with 1st array index = cell no, 2nd array being the genotype value
#print(transposed_df)

timepoint_cells, noOfCells = get_timepoint_cells(args.cells)
cell_genotype = transposed_df.values.tolist()
#print(" cell genotype ",cell_genotype)

# Loop over timepoint_cells and for each cell get the genotype and the mutation based on the index of the genotype array.
get_timepoint_mutations(timepoint_cells, cell_genotype, args.op_path)
