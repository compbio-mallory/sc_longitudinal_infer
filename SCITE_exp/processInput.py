import pandas as pd
import argparse

def process_input(input_file, mut_file, cell_file, output_file, text_output):
    input_df = pd.read_csv(input_file, sep='\t', header=None)
    transposed_df = input_df.T
    transposed_df.to_csv(output_file, sep=' ', index=False, header=False)

    mfile = open(mut_file,"r")
    mline = mfile.readline().rstrip('\n')
    mutSet = set()
    while(mline != ""):
        mline = mline.split("\t")
        mutations = mline[1].split(';')
        mutSet.update(mutations)
        mline = mfile.readline().rstrip('\n')
    mutCount = len(mutSet)
    print(" Total mutation count ",mutCount)

    cfile = open(cell_file,"r")
    cline = cfile.readline().rstrip('\n')
    cellCount = 0
    while(cline != ""):
        cline = cline.split("\t")
        cells = cline[1].split(';')
        cellCount = cellCount+len(cells)
        cline = cfile.readline().rstrip('\n')
    print(" Total cell count ",cellCount)

    scite_op_path_1 = input_file.split('/')[1]
    scite_op_path_2 = input_file.split('/')[2]
    scite_op_path_3 = input_file.split('/')[3]
    print(" SCITE op path ",scite_op_path_1," ",scite_op_path_2," ",scite_op_path_3)
    with open(text_output,"w+") as f:
        scite_str = "../../SCITE/scite -i "+output_file+" -n "+str(mutCount)+" -m "+str(cellCount)+" -r 1 -l 1565000 -fd 0.01 -ad 0.2 -a -max_treelist_size 1 -o ./sim_output_t3/"+scite_op_path_2+"/"+scite_op_path_3+"/"+scite_op_path_2+"_"+scite_op_path_3
        #print(scite_str)
        f.write(scite_str)
        f.write("\n")
        f.write("#")
        f.write("\n")

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input data matrix")
parser.add_argument("-mut", "--mut",dest ="mut", help="No. of mutations")
parser.add_argument("-cells", "--cells",dest ="cells", help="No. of cells")
parser.add_argument("-output", "--output",dest ="output", help="Output data matrix")
parser.add_argument("-txtOp", "--txtOp",dest ="txtOp", help="Output text file to save the command to run SCITE")
args = parser.parse_args()

process_input(args.input, args.mut, args.cells, args.output, args.txtOp)
