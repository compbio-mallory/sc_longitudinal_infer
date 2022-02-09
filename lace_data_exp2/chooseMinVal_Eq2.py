import cellToCluster
import glob
import os
import shutil
import subprocess
import argparse
import json

def call_equation2(inputFile, gp, cp, lambda_coeff, clusters, cellToCluster):
    eq2_cmd = 'python equation2Value.py -input '+inputFile+' -gp '+gp+' -cp '+cp+' -lambda_coeff '+lambda_coeff+' -subclones '+clusters+' -cellCluster '+cellToCluster
    print(eq2_cmd)
    result = subprocess.run(eq2_cmd, stdout=subprocess.PIPE, shell=True)
    return result.stdout.decode('utf-8')

def find_minEq2Value(eq2Val_dict):
   return min(eq2Val_dict, key=eq2Val_dict.get) 

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input matrix (similar to input to SCG)")
parser.add_argument("-lambda_coeff", "--lambda_coeff",dest = "lambda_coeff", help="Lambda co-efficient of Eq. 2")
args = parser.parse_args()

eq2Val_dict = {}
for dirname in glob.iglob('./scg_clusterNo*'):
    print(dirname)
    clusterNo = dirname.split('_')[2]
    cp = ""
    gp = ""
    #print(clusterNo)
    for filename in os.scandir(dirname):
        if filename.path.endswith('.tsv.gz'):
            #print(dirname+"/"+str(filename.name))
            if 'cluster' in filename.name:
                cp = dirname+"/"+filename.name
                cellToCluster.cell_cluster(dirname+"/"+str(filename.name),dirname+'/cellToCluster.json')
            else:
                gp = dirname+"/"+filename.name    
    #print(" CP "+cp+" GP "+gp)
    eq2_value = call_equation2(args.input, gp, cp, args.lambda_coeff, clusterNo, dirname+'/cellToCluster.json')
    eq2Val_dict[dirname] = eq2_value

with open('eq2_Values.json','w') as fp:
    json.dump(eq2Val_dict, fp)

minEq2_dir = find_minEq2Value(eq2Val_dict)
print(minEq2_dir) # Final output to select the cluster with min value of Eq 2.

# Remove other directories other than the minimum one
for dirname in glob.iglob('./scg_clusterNo*'):
    if minEq2_dir == dirname:
        continue
    else:
        subprocess.run('rm -r '+dirname, shell=True)

#print(eq2Val_dict)
