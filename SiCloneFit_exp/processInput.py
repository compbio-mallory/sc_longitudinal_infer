import pandas as pd
import argparse

def process_input(input_file, mut_file, cell_file, output_file, text_output):
    input_df = pd.read_csv(input_file, sep='\t', header=None)
    transposed_df = input_df.T
    transposed_df.to_csv(output_file, sep=' ', header=False)
    simDataType = input_file.split('/')[1]
    repk = input_file.split('/')[2]

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

    #op_fname = (text_output.split('/')[1]).split('.')[0]
    #cellnames_str = ['sc'+str(i) for i in range(1,cellCount+1)]
    #print("Cellnames ",cellnames_str)
    #with open('processed_sim_input/'+op_fname+'_cellNames.txt','w+') as f:
    #    f.write(' '.join(cellnames_str))
    #    f.write("\n")

    #with open('processed_sim_input/'+op_fname+'_geneNames.txt','w+') as f:
    #    for i in range(1,mutCount+1):
    #        f.write(str(i))
    #        f.write("\n")

    with open(text_output,'w+') as f:
        siclonefit_run = 'java -jar ../../siclonefit/SiCloneFiTComplete.jar -m '+str(cellCount)+' -n '+str(mutCount)+' -fp 0.01 -fn 0.2 -r 10 -iter 100 -df 0 -ipMat '+output_file+' -outDir sim_output_t3_trees/'+simDataType+'/'+repk
        #sifit_inferMut = 'java -cp ../../SiFit/SiFit.jar SiFit.algorithm.InferAncestralStates -fp 0.01 -fn 0.2 -df 0 -ipMat '+output_file+' -tree processed_sim_input/'+op_fname+'.D_mlTree.newick -cellNames processed_sim_input/'+op_fname+'_cellNames.txt -geneNames processed_sim_input/'+op_fname+'_geneNames.txt -expectedMatrix inferredGenotype/'+op_fname+'.csv'
        f.write(siclonefit_run)
        f.write("\n")
        #f.write(sifit_inferMut)
        #f.write("\n")
        f.write("#")
        f.write("\n")

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input data matrix")
parser.add_argument("-mut", "--mut",dest ="mut", help="No. of mutations")
parser.add_argument("-cells", "--cells",dest ="cells", help="No. of cells")
parser.add_argument("-output", "--output",dest ="output", help="Output data matrix")
parser.add_argument("-txtOp", "--txtOp",dest ="txtOp", help="Output text file to save the command to run SiFit")
args = parser.parse_args()

process_input(args.input, args.mut, args.cells, args.output, args.txtOp)
