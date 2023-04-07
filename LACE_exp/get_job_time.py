''' This script saves the time for each slurm job from the jobID. To plot the boxplots we need to read the Wall clock only for completed jobs. '''
import subprocess

''' Get the slurm job ids from HPC. '''
def getJobIDs():
    getJobIDcmd = "sacct -S 2023-03-29T10:00:00 > ./jobIds.txt"
    subprocess.call(getJobIDcmd, shell=True)
    jobID_file = open("jobIds.txt",'r')
    jf_line = jobID_file.readline().rstrip()
    jobIds_list = {}
    while(jf_line != ""):
        if "JobID" in jf_line or "bat+" in jf_line or "ext+" in jf_line or '------------' in jf_line or 'LACE' not in jf_line:
            jf_line = jobID_file.readline().rstrip()
            continue
        jf_line_arr = jf_line.split(" ")

        #if "COMPLETED" not in jf_line_arr or "CANCELLED" not in jf_line_arr: # We include Completed and Cancelled job info and handle the time in boxplot script
        #    jf_line = jobID_file.readline().rstrip()
        #    continue

        #print(jf_line_arr[0])
        jobID = int(jf_line_arr[0])
        if jf_line_arr[10] == '':
            jobName = jf_line_arr[11]
        else:
            jobName = jf_line_arr[10]
        jobIds_list[jobID] = jobName
        jf_line = jobID_file.readline().rstrip()
    return jobIds_list

# Save the time for each jobs.
def save_jobTime(jobIds_list):
    for jobid, jobname in jobIds_list.items():
        slurm_job_no = jobname.split('E')[1]
        slurm_file = "sim_run_LACE."+slurm_job_no+".slurm"
        print(" Job ID ",jobid," Slurm file ",slurm_file)
        sf = open(slurm_file,'r')
        sline = sf.readline().rstrip('\n')
        while(sline != ""):
            if "SBATCH" in sline and "--output" in sline:
                sline = ((sline.split(" "))[1]).split("/")
                #if "cov03" in sline or "lowCoverage" in sline:
                #    slurm_job_no = slurm_job_no+1
                #    continue
                op_dir = "sim_output/"+sline[1]+"/"+sline[2]+"/"
                #print(" Output dir ",op_dir)
                #f_op = open(op_dir+"/time.txt",'w')
                #print(f_op)
                job_stats = "seff "+str(jobid)+" > "+op_dir+"/time.txt"
                #subprocess.call([job_stats," >> ",op_dir+"/time.txt"], shell=True)
                subprocess.call(job_stats, shell=True)

            sline = sf.readline().rstrip('\n')

jobIds_list = getJobIDs()
print("No. of jobs ",jobIds_list)
save_jobTime(jobIds_list)


