import json
import glob
import os

for dirname in glob.iglob('./scg_clusterNo*'):
    print(dirname)
    #print(clusterNo)
    clusterToCell_dict = {}
    for filename in os.scandir(dirname):
        if filename.path.endswith('.json'):
            print(filename.name)
            with open(dirname+"/"+filename.name,'r') as fp:
                cellToCluster_dict = json.loads(fp.read())
                #print(cellToCluster_dict)
                for key, value in cellToCluster_dict.items():
                    if value in clusterToCell_dict:
                        tmp_list = clusterToCell_dict[value]
                        tmp_list.append(key)
                        clusterToCell_dict[value] = tmp_list
                    else:
                        cell_list = []
                        cell_list.append(key)
                        clusterToCell_dict[value] = cell_list
    print(clusterToCell_dict)
    with open('clusterToCell.json','w') as fp:
        json.dump(clusterToCell_dict,fp,indent=4)
