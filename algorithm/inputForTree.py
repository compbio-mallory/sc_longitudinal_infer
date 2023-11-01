import argparse
import numpy as np

''' For each cluster in the assignment file get the cells. '''
def get_cell_cluster(assignment):
    with open(assignment,'r') as f:
        lines = f.readlines()

    assignment_list = ((lines[1].split('\t'))[2]).split(" ")
    cell_cluster = {} # key is clusterID, value is cell ID.
    #print(len(assignment_list))
    cell_count = 0
    for i in assignment_list:
        if '\n' in i:
            i = i.replace("\n","")
        if i in cell_cluster:
            temp_cell_list = cell_cluster[i]
            temp_cell_list.append(cell_count)
            cell_cluster[i] = temp_cell_list
            cell_count = cell_count+1
        else:
            temp_cell_list = []
            temp_cell_list.append(cell_count)
            cell_cluster[i] = temp_cell_list
            cell_count = cell_count+1

    return cell_cluster

''' Get the cells of each cluster from the simulated file. '''
def sim_cell_cluster(clusterFile):
    with open(clusterFile,'r') as f:
        lines = f.readlines()

    cell_cluster = {}
    for line in lines:
        #print(line.split('\t'))
        clusterID = ((line.split('\t')[0]).split('_'))[1]
        cells = (line.split('\t')[1].strip()).split(';')
        cell_cluster[clusterID] = cells

    return cell_cluster

''' Get timepoints for each cell from the simulation file. '''
def get_cell_timepoints(celltp_file):
    file = open(celltp_file,"r")
    line = file.readline().rstrip('\n')
    cell_tp = {}
    while(line != ""):
        line_a = line.split('\t')
        #print(line_a)
        tp = (line_a[0].split('_'))[0]
        cells = line_a[1].split(';')
        for cell in cells:
            cell_tp[cell] = tp
        line = file.readline().rstrip('\n')
    return cell_tp

''' Get cluster genotype by taking a majority vote for each cell. '''
def vote_cluster_genotypes(cell_gs_list):
    cell_gs_len = len(cell_gs_list)
    cg_len = len(cell_gs_list[0])
    #print(cell_gs_len," ",cg_len)

    #print(cell_gs_list[0])
    #print("============================")
    voted_cg = []

    for j in range(0,cg_len):
        total0s = 0
        total1s = 0
        for i in range(0,cell_gs_len):
            if cell_gs_list[i][j] == '0':
                total0s = total0s+1
            if cell_gs_list[i][j] == '1':
                total1s = total1s+1
        #print(" Total 0s ",total0s," Total 1s ",total1s)
        if total0s > total1s:
            voted_cg.append('0')
        if total1s >= total0s:
            voted_cg.append('1')
    return voted_cg

''' Get the cluster genotype or G' matrix. '''
def get_cluster_genotype(cell_cluster, consensus_genotype):
    cluster_genotype = {}
    cell_no = 0
    firstLine = True
    with open(consensus_genotype) as file:
        for line in file:
            if firstLine == True:
                firstLine = False
                continue
            for cluster, cell in cell_cluster.items():
                cell_val = (cell_cluster[cluster])
                for cell in cell_val:
                    if cell_no == cell:
                        print(cluster," --cell ",cell)
                        l = (line.strip()).split('\t')
                        print(l[1:])
                        if cluster in cluster_genotype:
                            temp_list = cluster_genotype[cluster]
                            temp_list.append(l[1:])
                            cluster_genotype[cluster] = temp_list
                        else:
                            temp_list = []
                            temp_list.append(l[1:])
                            cluster_genotype[cluster] = temp_list
            cell_no = cell_no+1
    #print(cluster_genotype.get('2'))
    for cluster, cell_gs in cluster_genotype.items():
        #print(" Cluster ",cluster)
        if len(cell_gs) == 1:
            cluster_genotype[cluster] = cell_gs[0]
            continue
        voted_cg = vote_cluster_genotypes(cell_gs) 
        cluster_genotype[cluster] = voted_cg
    return cluster_genotype

''' Get the G' matrix from simulated ground truth. '''
def sim_cluster_genotype(cell_cluster, genotype):
    cluster_genotype = {}
    cell_no = 0
    with open(genotype) as file:
        for line in file:
            for cluster, cell in cell_cluster.items():
                cell_val = (cell_cluster[cluster])
                for cell in cell_val:
                    if cell_no == int(cell):
                        l = (line.strip()).split('\t')
                        if cluster in cluster_genotype:
                            temp_list = cluster_genotype[cluster]
                            temp_list.append(l)
                            cluster_genotype[cluster] = temp_list
                        else:
                            temp_list = []
                            temp_list.append(l)
                            cluster_genotype[cluster] = temp_list
            cell_no = cell_no+1

    for cluster, cell_gs in cluster_genotype.items():
        cluster_genotype[cluster] = cell_gs[0]
        #if len(cell_gs) == 1:
        #    cluster_genotype[cluster] = cell_gs[0]
        #    continue
        #voted_cg = vote_cluster_genotypes(cell_gs)
        #cluster_genotype[cluster] = voted_cg
    return cluster_genotype

''' Merge clusters with same genotypes. Update the cluster_genotype and cell_cluster. '''
def merge_clusters(cluster_genotype, cell_cluster):
    same_clusters = {} # genotype will be the key and clusters values.
    for cluster, genotype in cluster_genotype.items():
        for c1, g1 in cluster_genotype.items():
            if cluster == c1:
                continue
            if genotype == g1:
                print(" Same ",g1)
                key_g1 = ",".join(g1)
                if key_g1 in same_clusters:
                    temp_list = same_clusters[key_g1]
                    temp_list.append(c1)
                    temp_list.append(cluster)
                    same_clusters[key_g1] = temp_list
                else:
                    temp_list = []
                    temp_list.append(c1)
                    temp_list.append(cluster)
                    same_clusters[key_g1] = temp_list
    print(" Same clusters ",same_clusters)

    for genotype,cluster in same_clusters.items():
        set_cluster = set(cluster)
        merged_cells = []
        for c in set_cluster:
            cells = cell_cluster[c]
            merged_cells.extend(cells)
            del cell_cluster[c]
            del cluster_genotype[c]
            keep_cluster = c
        print(" Merged cells ",merged_cells)
        cell_cluster[keep_cluster] = merged_cells
        cluster_genotype[keep_cluster] = genotype.split(",")
    return cluster_genotype, cell_cluster

def saveGmatrix(cell_cluster, cluster_genotype, cell_tp, out_file):
    out_f = open(out_file, 'w')
    for cluster,cells in cell_cluster.items():
        timepoint = set()
        for cell in cells:
            #print(cell)
            #print(cell_tp)
            timepoint.add(cell_tp[str(cell)])
        for tp in timepoint:
            key = cluster+"_"+tp
            print(cluster_genotype[cluster])
            value = "\t".join(cluster_genotype[cluster])
            print("\t".join([str(key),str(value)]), file = out_f)
    out_f.close()

''' Save simulated G matrix here. '''
def saveSimGmatrix(cell_cluster, cluster_genotype, cell_tp, out_file):
    out_f = open(out_file, 'w')
    for cluster,cells in cell_cluster.items():
        timepoint = set()
        for cell in cells:
            #print(cell)
            #print(cell_tp)
            timepoint.add(cell_tp[cell])
        for tp in timepoint:
            key = cluster+"_"+tp
            print(cluster_genotype[cluster])
            value = "\t".join(cluster_genotype[cluster])
            print("\t".join([str(key),str(value)]), file = out_f)
    out_f.close()

parser = argparse.ArgumentParser()
parser.add_argument("-cc", "--cc",dest = "cc", help="Assignment.txt file indicating clusters of cells.")
parser.add_argument("-cg", "--cg",dest = "cg", help="Consensus genotype from BnpC.")
parser.add_argument("-celltp", "--celltp",dest = "celltp", help="Cell timepoints.")
parser.add_argument("-sim", "--sim",dest = "sim", help="Simulated cell assignment file.")
parser.add_argument("-op", "--op",dest = "op", help="Output file to save G prime matrix.")
parser.add_argument("-ccop","--ccop",dest = "ccop", help="Output file to save the cell cluster dictionary.")
args = parser.parse_args()

if args.sim == "true":
    cell_cluster = sim_cell_cluster(args.cc)
else:
    cell_cluster = get_cell_cluster(args.cc)
print(" Cell cluster ==== ")
print(cell_cluster)

if args.sim == "true":
    cluster_genotype = sim_cluster_genotype(cell_cluster, args.cg)
else:
    cluster_genotype = get_cluster_genotype(cell_cluster, args.cg)

print(len(cluster_genotype))
cell_tp = get_cell_timepoints(args.celltp)
#print(" Cell timepoints ",cell_tp)
if args.sim == "true":
    saveSimGmatrix(cell_cluster,cluster_genotype,cell_tp,args.op)
else:
    cluster_genotype, cell_cluster = merge_clusters(cluster_genotype, cell_cluster)
    print(" UPDATED ====== ")
    print(cluster_genotype)
    print(cell_cluster)
    np.save(args.ccop, cell_cluster) # Save the cell cluster file and read it later.
    # Load
    #read_dictionary = np.load('my_file.npy',allow_pickle='TRUE').item()
    #print(read_dictionary['0'])
    saveGmatrix(cell_cluster,cluster_genotype,cell_tp,args.op)
