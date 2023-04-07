import pandas as pd
import argparse

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input D matrix")
parser.add_argument("-output", "--output",dest ="output", help="Output to save the results")
args = parser.parse_args()

input_df = pd.read_csv(args.input, sep='\t', header=None)
transposed_df = input_df.T

#print(transposed_df)
print(" Output path ",args.output)
#transposed_df.to_csv("real_data_input/real_data_for_bnpc.tsv", sep='\t', index=False, header=False)

# Save input for BnpC 
op_fname = (((args.input.split('/'))[4]).split('.'))[0]
print(op_fname)
transposed_df.to_csv(args.output+'/'+op_fname+'.D.tsv', sep='\t', index=False, header=False)
