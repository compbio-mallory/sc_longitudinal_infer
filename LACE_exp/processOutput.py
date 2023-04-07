import argparse
import json

def getCellsInClones(cInf_json,assignment_path):
    f = open(cInf_json,'r')
    cInf_data = json.load(f)
    clone_tp = {} # This dict saves the clones at each timepoint.
    for k,timepoints in cInf_data.items():
        tp_count = 1
        wholeTree_clone_cells = []
        for tp,clones in timepoints.items():
            clone_set = set()
            clone_cells = []
            for clone in clones:
                clone = clone[0]
                clone_cells.append(clone)
                wholeTree_clone_cells.append(clone)
                #print(" Clone ",clone[0])
                if tp_count in clone_tp:
                    temp_clone = clone_tp[tp_count]
                    temp_clone.add(clone)
                    clone_tp[tp_count] = temp_clone
                else:
                    clone_set.add(clone)
                    clone_tp[tp_count] = clone_set

            w_file = open(assignment_path+'/assignment_tp'+str(tp_count)+'.txt',"w") # assignment file to calculate the V-meas.
            w_file.write(' '.join(list(map(str,clone_cells))))
            tp_count = tp_count+1

        assignment_file = open(assignment_path+'/assignment.txt',"w")
        assignment_file.write(' '.join(list(map(str,wholeTree_clone_cells))))
    print(" Clones in timepoints ",clone_tp)
    print(" Whole tree clones ",len(wholeTree_clone_cells))
    return clone_tp

def getCloneSNVs(clonesSum_json,clone_tp,assignment_path):
    f = open(clonesSum_json,'r')
    cloneSum_data = json.load(f)
    clone_SNVs = {}
    tp_SNVs = {}
    # Gets the SNVs assigned to each clones
    for k,clones in cloneSum_data.items():
        for cl,snv in clones.items():
            snv_list = []
            if type(snv) == str:
                print(snv.split('V'))
                mut_no = snv.split('V')[1]
                snv_list.append(mut_no)
            else:
                for s in snv:
                    mut_no = s.split('V')[1]
                    snv_list.append(mut_no)
            clone_SNVs[cl] = snv_list

    # Get the mutations in each timepoint
    for clone,SNVs in clone_SNVs.items():
        mut_list = []
        for tp,clones in clone_tp.items():
            key = str(tp)+"_"+str(tp+1)
            for c in clones:
                cloneNo = clone.split('_')[1]
                if c == int(cloneNo):
                    if key in tp_SNVs:
                        temp_list = tp_SNVs[key]
                        temp_list.extend(SNVs)
                        tp_SNVs[key] = temp_list
                    else:
                        temp_list = []
                        temp_list.extend(SNVs)
                        tp_SNVs[key] = temp_list

    # Remove duplicates from the timepoint SNVs dictionary.
    for tp,snvs in tp_SNVs.items():
        snvs = set(snvs)
        tp_SNVs[tp] = list(snvs)

    w_file = open(assignment_path+'/tp_mut.json',"w")
    json.dump(tp_SNVs,w_file)
    print(" Mutations at each timepoint ",tp_SNVs)
    print("----------------------------")
    print(clone_SNVs)

parser = argparse.ArgumentParser()
parser.add_argument("-cInf", "--cInf",dest ="cInf", help="C inference results from LACE")
parser.add_argument("-clonesSum", "--clonesSum",dest ="clonesSum", help="Summary of mutations present in each clone of LACE")
parser.add_argument("-path", "--path",dest ="path", help="Path to save the assignment file")
args = parser.parse_args()

clone_tp = getCellsInClones(args.cInf,args.path)
getCloneSNVs(args.clonesSum,clone_tp,args.path)
