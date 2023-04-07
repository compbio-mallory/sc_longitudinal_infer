import argparse 

# Create job file for RobustClone.
def create_job_file(fname,n,mpc,p,cluster_algo,sim):
    filename = open(fname,'r')
    line = filename.readline().rstrip('\n')
    op_fname = (fname.split("."))[0]
    count = 0
    lines_together = []
    temp_line = []
    while(line != ""):
        if "#" in line:
            lines_together.append(temp_line)
            temp_line = []
        temp_line.append(line)
        line = filename.readline().rstrip('\n')
    #print(lines_together)
    for lt in lines_together:
        #print(lt)
        line_items = lt[len(lt)-1].split(" ")
        for item in line_items:
            if 'sim_output' in item:
                op_dir_list = item.split('/')
                op_dir = op_dir_list[1]+'/'+op_dir_list[2]
                print(" Output dir ",op_dir)
        with open(op_fname+"."+str(count)+".slurm",'w') as f:
            f.write("#!/bin/bash\n")
            f.write("#SBATCH --job-name="+cluster_algo+str(count)+"\n")
            f.write("#SBATCH -N 1\n")
            f.write("#SBATCH -n "+str(n)+"\n")
            f.write("#SBATCH --mem-per-cpu="+mpc+"\n")
            f.write("#SBATCH -p "+p+"\n")
            f.write("#SBATCH --mail-type=ALL\n")
            f.write("#SBATCH --output=sim_output/"+op_dir+"/output"+str(count)+".out\n")
            f.write("#SBATCH --error=sim_output/"+op_dir+"/error"+str(count)+".out\n")
            for st in lt:
                if st == "#":
                    continue
                f.write(st)
                f.write('\n')
        count=count+1

parser = argparse.ArgumentParser()
parser.add_argument("-input", "--input",dest ="input", help="Input .sh file")
parser.add_argument("-n","--n",dest="n", help="No of cpus")
parser.add_argument("-mem_per_cpu", "--mem_per_cpu",dest ="mem_per_cpu", help="Memory per CPU")
parser.add_argument("-p", "--p",dest ="p", help="partition")
parser.add_argument("-algo", "--algo",dest ="algo", help="Which cluster algo you are running")
parser.add_argument("-sim", "--sim",dest ="sim", help="Simulated data")
args = parser.parse_args()

create_job_file(args.input, args.n, args.mem_per_cpu, args.p, args.algo, args.sim)

