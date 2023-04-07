""" Read the file merged.tsv and create a file similar to SNVcell.csv. It represents the following: 
    timepoint_timepoint cell1;cell2;cell3 
    We can save the timepoints of the cells here in a csv file to be read by the script inputForTree.py """

def save_csv(input_file, output_file):
    ifile = open(input_file, 'r')
    ifile_line = ifile.readline().rstrip('\n')
    count = 0
    cell_count = 0
    cell_tp = {}
    while(ifile_line != ""):
        if count == 0:
            count = count+1
            ifile_line = ifile.readline().rstrip('\n')
        cell = ifile_line.split('\t')[0]
        if "X1" in cell:
            timepoint = "2_2"
        elif "X2" in cell:
            timepoint = "3_3"
        elif "X4" in cell:
            timepoint = "4_4"

        if timepoint in cell_tp:
            cell_list = cell_tp[timepoint]
            cell_list.append(str(cell_count))
            cell_tp[timepoint] = cell_list
            cell_count = cell_count+1
        else:
            cell_list = []
            cell_list.append(str(cell_count))
            cell_tp[timepoint] = cell_list
            cell_count = cell_count+1
        ifile_line = ifile.readline().rstrip('\n')
    print(" Cell timepoints ",cell_tp)

    with open(output_file, 'w') as f:
        for k,v in cell_tp.items():
            f.write(k+'\t'+';'.join(v))
            f.write('\n')

save_csv("merged.tsv","cell_timepoints.csv")


