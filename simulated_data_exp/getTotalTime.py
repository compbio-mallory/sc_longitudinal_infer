import argparse
import os
import math

def totalTime(inputDir):
    print(" Dir name ",inputDir)
    total_time = 0
    for subdir, dirs, files in os.walk(inputDir):
        print(" subdir ",subdir)
        print(" files ",files)
        subdir_len = len(subdir.split("/"))
        for f in files:
            if "lgError" in f: # Get 
                eFilePath = os.path.join(subdir, f)
                time_object = 0
                with open(eFilePath,'r') as f1:
                    lines = f1.readlines()
                    #print(lines)
                    for lgFile_line in lines:
                        if "real" in lgFile_line:
                            time_split = (((lgFile_line.split('\t'))[1]).strip()).split('m')
                            #print(" Time split ",time_split)
                            minute_time = int(time_split[0])
                            if minute_time != 0:
                                minute_time = minute_time*60
                            second_time = int(math.ceil(float((time_split[1].split('s'))[0])))
                            time_object = time_object+minute_time+second_time
                            print(" Time objects ",time_object)

                total_time = total_time+time_object
                #print(" Total time ",total_time)

        if subdir_len == 4: # Get the time of running BnpC
            for f in files:
                if "error" in f and "out" in f:
                    eFilePath = os.path.join(subdir, f)
                    with open(eFilePath,'r') as f1:
                        lines = f1.readlines()
                        for bnpcLine in lines:
                            if "real" in bnpcLine:
                                time_split = (((bnpcLine.split('\t'))[1]).strip()).split('m')
                                #print(" Time split ",time_split)
                                minute_time = int(time_split[0])
                                #print(" Mins time ",minute_time)
                                if minute_time != 0:
                                    minute_time = minute_time*60
                                second_time = int(math.ceil(float((time_split[1].split('s'))[0])))
                                #print(" Total time here ",total_time)
                                #print(" Before summation time ",minute_time+second_time)
                                total_time = total_time+minute_time+second_time
                                #print(" Total time after sum ",total_time)

        #total_time = total_time+time_object
        print(" Total time ",total_time)

    with open(inputDir+'/'+'total_time.txt','w') as of: # Save the total time in a .txt file to be read for boxplots
        of.write(str(total_time))
        of.write('\n')

parser = argparse.ArgumentParser()
parser.add_argument("-sim_dir","--sim_dir",dest="sim_dir", help="The directory where the files with running time are present.")
args = parser.parse_args()
totalTime(args.sim_dir)
