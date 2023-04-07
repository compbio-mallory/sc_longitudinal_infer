import argparse
import pandas as pd

def read_file(input_file, output_file):
    rfile = open(input_file,'r')
    rfile_line = rfile.readline().rstrip('\n')
    count = 0
    mut_list = []
    while(rfile_line != ""):
        if count == 0:
            count = count+1
            rfile_line = rfile.readline().rstrip('\n')
            continue
        mut_line = rfile_line.split('\t')[1:]
        mut_list.append(mut_line)
        rfile_line = rfile.readline().rstrip('\n')
    print("No of mutations ",len(mut_list[0]))
    #print(" Mutation lists ",mut_list)
    df = pd.DataFrame(mut_list)
    df.to_csv('input.D.csv',sep='\t', index=False, header=False)
    df = df.T
    print(" Pandas dataframe ",df)
    df.to_csv(output_file+'.D.tsv', sep='\t', index=False, header=False)

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input real data matrix")
parser.add_argument("-output", "--output", dest="output", help="Output file to save the data")
args = parser.parse_args()
read_file(args.input, args.output)

