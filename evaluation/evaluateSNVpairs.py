import ast
import argparse

def read_files(SNVfile):
    # check if the saved file can be read
    print(SNVfile)
    with open(SNVfile,'r') as f:
        #print(f.readline().rstrip('\n'))
        #if f.readline().rstrip('\n') == set():
        #    print(True)
        #    SNV_set = set()
        #else:
            #print(False)
        try:
            SNV_set = ast.literal_eval(f.read())
        except:
            SNV_set = set()
    return SNV_set

''' Convert sets to lists to check SNV pairs '''
def convertSetsToLists(givenSet):
    convertedList = []
    for elem in givenSet:
        elem_list = elem.split('_')
        convertedList.append(elem_list)
    #print(convertedList)
    return convertedList

''' Correct number of SNVs when compared with ground truth. '''
#def correctSNVs(gt_list, other_list, SNVtype):
#    correct_pairs = 0
#    gt_list_len = len(gt_list)
#    for g in gt_list:
#        for o in other_list:
#            if (g[0] == o[0] or g[0] == o[1]) and (g[1] == o[1] or g[1] == o[0]):
            #if g[0] == o[0] and g[1] == o[1]:
#                     correct_pairs = correct_pairs+1
                     #print(g," ---- ",o)
#    print(" Correct pairs ",correct_pairs," GT len ",gt_list_len)
#    if gt_list_len == 0:
#        correctSNVs_perc = 0
#    else:
#        correctSNVs_perc = (correct_pairs/gt_list_len) * 100
#    print(SNVtype+" = "+str(round(correctSNVs_perc)))
#    return correct_pairs, gt_list_len

''' Correct number of SNVs when compared with ground truth. '''
def correctSNVs(gt_list, other_list, SNVtype):
    correct_pairs = 0
    gt_list_len = len(gt_list)
    for g in gt_list:
        possible_pair1 = [g[0]]
        possible_pair1.append(g[1])
        possible_pair2 = [g[1]]
        possible_pair2.append(g[0])
        if possible_pair1 in other_list or possible_pair2 in other_list:
            correct_pairs = correct_pairs+1
    print(" Correct pairs ",correct_pairs," GT len ",gt_list_len)
    if gt_list_len == 0:
        correctSNVs_perc = 0
    else:
        correctSNVs_perc = (correct_pairs/gt_list_len) * 100
    print(SNVtype+" = "+str(round(correctSNVs_perc)))
    return correct_pairs, gt_list_len
#def compareSameSNVs(gt_sameSNVs_list, other_sameSNVs_list):

#def compareAncestralSNVs(gt_ancestralSNVs_list, other_ancestralSNVs):

#def compareParallelSNVs(gt_parallelSNVs, other_parallelSNVs):

parser = argparse.ArgumentParser()
parser.add_argument("-gt", "--gt",dest ="gt", help="Ground truth file path for SNVs")
parser.add_argument("-siclonefit", "--siclonefit",dest ="siclonefit", help="SiCloneFit file path for SNVs")
parser.add_argument("-scLongTree", "--scLongTree",dest ="scLongTree", help="scLongTree file path for SNVs")
parser.add_argument("-scite", "--scite",dest ="scite", help="SCITE file path for SNVs")
parser.add_argument("-lace", "--lace",dest ="lace", help="LACE file path for SNVs")
args = parser.parse_args()

# Get all the SNV files
gt_sameSNVs = read_files(args.gt+"/sameSNVs.txt")
gt_ancestralSNVs = read_files(args.gt+"/ancestralSNVs.txt")
gt_parallelSNVs = read_files(args.gt+"/inComSNVs.txt")

try:
    siclonefit_sameSNVs = read_files(args.siclonefit+"/sameSNVs.txt")
    siclonefit_ancestralSNVs = read_files(args.siclonefit+"/ancestralSNVs.txt")
    siclonefit_parallelSNVs = read_files(args.siclonefit+"/inComSNVs.txt")
except:
        print(" SiCloneFit file not found. ")

scLongTree_sameSNVs = read_files(args.scLongTree+"/sameSNVs.txt")
scLongTree_ancestralSNVs = read_files(args.scLongTree+"/ancestralSNVs.txt")
scLongTree_parallelSNVs = read_files(args.scLongTree+"/inComSNVs.txt")

try:
    scite_sameSNVs = read_files(args.scite+"/sameSNVs.txt")
    scite_ancestralSNVs = read_files(args.scite+"/ancestralSNVs.txt")
    scite_parallelSNVs = read_files(args.scite+"/inComSNVs.txt")
except:
    print(" SCITE file doesn't exist. ")

try:
    lace_sameSNVs = read_files(args.lace+"/sameSNVs.txt")
    lace_ancestralSNVs = read_files(args.lace+"/ancestralSNVs.txt")
    lace_parallelSNVs = read_files(args.lace+"/inComSNVs.txt")
except:
        print(" LACE file not found. ")

#print(gt_parallelSNVs)
gt_sameSNVs_lst = convertSetsToLists(gt_sameSNVs)
gt_parallelSNVs_lst = convertSetsToLists(gt_parallelSNVs)
gt_ancestralSNVs_lst = convertSetsToLists(gt_ancestralSNVs)

try:
    siclonefit_sameSNVs_lst = convertSetsToLists(siclonefit_sameSNVs)
    siclonefit_ancestralSNVs_lst = convertSetsToLists(siclonefit_ancestralSNVs)
    siclonefit_parallelSNVs_lst = convertSetsToLists(siclonefit_parallelSNVs)
except:
        print(" SiCloneFit doesn't exist for this. ")

scLongTree_sameSNVs_lst = convertSetsToLists(scLongTree_sameSNVs)
scLongTree_ancestralSNVs_lst = convertSetsToLists(scLongTree_ancestralSNVs)
scLongTree_parallelSNVs_lst = convertSetsToLists(scLongTree_parallelSNVs)

try:
    scite_sameSNVs_lst = convertSetsToLists(scite_sameSNVs)
    scite_ancestralSNVs_lst = convertSetsToLists(scite_ancestralSNVs)
    scite_parallelSNVs_lst = convertSetsToLists(scite_parallelSNVs)
except:
    print(" SCITE file doesn't exist. ")

try:
    lace_sameSNVs_lst = convertSetsToLists(lace_sameSNVs)
    lace_ancestralSNVs_lst = convertSetsToLists(lace_ancestralSNVs)
    lace_parallelSNVs_lst = convertSetsToLists(lace_parallelSNVs)
except:
        print(" LACE doesn't exist for this.")


print(" SiCloneFit's evaluation -------------- ")
try:
    correct_sameSNVs, gt_sameSNVs = correctSNVs(gt_sameSNVs_lst, siclonefit_sameSNVs_lst, "sameSNVs")
    correct_ancesSNVs, gt_ancesSNVs = correctSNVs(gt_ancestralSNVs_lst, siclonefit_ancestralSNVs_lst, "ancestralSNVs")
    correct_inComSNVs, gt_inComSNVs = correctSNVs(gt_parallelSNVs_lst, siclonefit_parallelSNVs_lst, "inComparableSNVs")
    print("siCloneFit SNV eval", (correct_sameSNVs+correct_ancesSNVs+correct_inComSNVs)/(gt_sameSNVs+gt_ancesSNVs+gt_inComSNVs))
except:
        print(" No SiCloneFit's evaluation")

print(" ScLongTree's evaluation -------------- ")
correct_sameSNVs, gt_sameSNVs = correctSNVs(gt_sameSNVs_lst, scLongTree_sameSNVs_lst, "sameSNVs")
#print(" scLongTree's ancestralSNVs ")
correct_ancesSNVs, gt_ancesSNVs = correctSNVs(gt_ancestralSNVs_lst, scLongTree_ancestralSNVs_lst, "ancestralSNVs")
correct_inComSNVs, gt_inComSNVs = correctSNVs(gt_parallelSNVs_lst, scLongTree_parallelSNVs_lst, "inComparableSNVs")
print("scLongTree SNV eval", (correct_sameSNVs+correct_ancesSNVs+correct_inComSNVs)/(gt_sameSNVs+gt_ancesSNVs+gt_inComSNVs))

print(" SCITE's evaluation -------------- ")
try:
    correct_sameSNVs, gt_sameSNVs = correctSNVs(gt_sameSNVs_lst, scite_sameSNVs_lst, "sameSNVs")
    correct_ancesSNVs, gt_ancesSNVs = correctSNVs(gt_ancestralSNVs_lst, scite_ancestralSNVs_lst, "ancestralSNVs")
    correct_inComSNVs, gt_inComSNVs = correctSNVs(gt_parallelSNVs_lst, scite_parallelSNVs_lst, "inComparableSNVs")
    print("SCITE SNV eval", (correct_sameSNVs+correct_ancesSNVs+correct_inComSNVs)/(gt_sameSNVs+gt_ancesSNVs+gt_inComSNVs))
except:
    print(" No SCITE's evaluation ")

print(" LACE's evaluation -------------- ")
try:
    correct_sameSNVs, gt_sameSNVs = correctSNVs(gt_sameSNVs_lst, lace_sameSNVs_lst, "sameSNVs")
    correct_ancesSNVs, gt_ancesSNVs = correctSNVs(gt_ancestralSNVs_lst, lace_ancestralSNVs_lst, "ancestralSNVs")
    correct_inComSNVs, gt_inComSNVs = correctSNVs(gt_parallelSNVs_lst, lace_parallelSNVs_lst, "inComparableSNVs")
    print("LACE SNV eval", (correct_sameSNVs+correct_ancesSNVs+correct_inComSNVs)/(gt_sameSNVs+gt_ancesSNVs+gt_inComSNVs))
except:
        print(" No LACE's evaluation")

