import argparse
import numpy as np

def readMatrix(matrix_file, G_file):
    mfile = open(matrix_file,'r')
    mfile_line = mfile.readline().rstrip('\n')
    matrix_list = []
    count=0
    while(mfile_line != ""):
        if G_file == True:
            if count == 0:
                mfile_line = mfile.readline().rstrip('\n')
                count = count+1
                continue
            line_list = mfile_line.split('\t')[1:]
            matrix_list.append(line_list)
            mfile_line = mfile.readline().rstrip('\n')
        else:
            line_list = mfile_line.split('\t')
            matrix_list.append(line_list)
            mfile_line = mfile.readline().rstrip('\n')

    #print(" Matrix ",len(matrix_list[0]))
    #print(" ------------------ ")
    #print(" Matrix ",len(matrix_list[1]))
    return matrix_list

def calc_FPFN(G_list, D_list, op_file):
    est_FP = 0
    est_FN = 0
    TN = 0
    TP = 0
    for i in range(len(G_list)):
        for j in range(len(G_list[0])):
            if G_list[i][j] == '1' and D_list[i][j] == '0': # G matrix here is the inferred genotype from BnpC and we are comparing it with D.
                est_FN = est_FN+1
            if G_list[i][j] == '0' and D_list[i][j] == '1':
                est_FP = est_FP+1
            if D_list[i][j] == '0' and G_list[i][j] == '0':
                TN = TN+1
            if D_list[i][j] == '1' and G_list[i][j] == '1':
                TP = TP+1

    print(" Total estimated FP ",est_FP," True negatives ",TN)
    print(" Total estimated FN ",est_FN," True positives ",TP)
    #print(" FP rate deno ",est_FP+TP)
    #print(" FN rate deno ",est_FN+TN)

    est_FP_rate = est_FP/(est_FP+TP)
    est_FN_rate = est_FN/(est_FN+TN)
    rounded_est_FP_rate = round(est_FP_rate,5)
    rounded_est_FN_rate = round(est_FN_rate,2)
    print(" Estimated FP rate ",rounded_est_FP_rate)
    print(" Estimated FN rate ",rounded_est_FN_rate)
    # Save the FP FN values in a dictionary and save that as a numpy object
    FP_FN_dict = {}
    FP_FN_dict['FP'] = rounded_est_FP_rate
    FP_FN_dict['FN'] = rounded_est_FN_rate
    np.save(op_file,FP_FN_dict)

parser = argparse.ArgumentParser()
parser.add_argument("-G","--G",dest = "G", help="The G matrix.")
parser.add_argument("-D","--D",dest = "D", help="The D matrix.")
parser.add_argument("-op","--op",dest = "op", help="The output path to save results.")
args = parser.parse_args()

G_list = readMatrix(args.G, True)
D_list = readMatrix(args.D, False)

#print(" G matrix ----- ",G_list)
#print(" D matrix ----- ",D_list)
calc_FPFN(G_list, D_list, args.op)
