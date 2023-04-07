import argparse

def get_header(mutations):
    header_arr = []
    header_arr.append('""')
    for i in range(len(mutations)):
        header_arr.append('"'+"SNV"+str(i)+'"')
    print(header_arr)
    return header_arr

def process_input(input_file, output_file):
    file = open(input_file,"r")
    line = file.readline().rstrip('\n')
    noOfMutations = line.split('\t')
    print(" No. of mutations ",len(noOfMutations))
    modified_lines = []
    header_arr = get_header(noOfMutations)
    header_arr_str = "\t".join(str(x) for x in header_arr)
    print(" Header ",header_arr_str)
    count=1 # Uncomment this line to process real data.
    while(line != ""):
        #if count == 1:
        #    line = file.readline().rstrip('\n')
        #    count=count+1
        #    continue
        #line = line.replace("\t",",")
        #print(" Line ",line)
        modified_lines.append(line)
        line = file.readline().rstrip('\n')
    with open(output_file,"w+") as f:
        f.write(header_arr_str)
        f.write("\n")
        #print(" No of cells ",len(modified_lines))
        for l in range(1,len(modified_lines)+1):
            f.write('"'+"cell"+str(l)+'"'+'\t'+modified_lines[l-1])
            #f.write(modified_lines[l-1])
            f.write("\n")

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input data matrix")
parser.add_argument("-output", "--output",dest ="output", help="Output data matrix")
#parser.add_argument("-mut","--mut",dest="mut", help="No of mutations")
args = parser.parse_args()

process_input(args.input, args.output)
