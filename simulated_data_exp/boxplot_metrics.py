import os
import argparse
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import math

def set_box_color(bp, color):
    plt.setp(bp['boxes'], color=color)
    plt.setp(bp['whiskers'], color=color)
    plt.setp(bp['caps'], color=color)
    plt.setp(bp['medians'], color=color)

def plot_multipleBp(algo_metric_dict, data_list, y_max, y_min, output):
    plt.style.use('seaborn')
    fig = plt.figure(figsize=(12,12))
    axes = fig.subplots(nrows=2, ncols=2)
    
    #axes.tick_params(axis="y",labelsize='small')
    #fig.suptitle('Beta splitting', fontsize=12) # Uncomment this line if you want the title

    #boxplot_ForMultipleMetric('Accuracy',algo_metric_dict.get('scg_acc'),algo_metric_dict.get('scclone_acc'),algo_metric_dict.get('bnpc_acc'),data_list,output,axes[0,0])
    boxplot_ForMultipleMetric('V-measure',algo_metric_dict.get('lgTree_vmeas'),algo_metric_dict.get('LACE_vmeas'),data_list,y_max,y_min,output,axes[0,0])
    boxplot_ForMultipleMetric('Time',algo_metric_dict.get('lgTree_time'),algo_metric_dict.get('LACE_time'),data_list,y_max,y_min,output,axes[0,1])
    boxplot_ForMultipleMetric('Precision',algo_metric_dict.get('lgTree_prec'),algo_metric_dict.get('LACE_prec'),data_list,y_max,y_min,output,axes[1,0])
    boxplot_ForMultipleMetric('Recall',algo_metric_dict.get('lgTree_recall'),algo_metric_dict.get('LACE_recall'),data_list,y_max,y_min,output,axes[1,1])
    #boxplot_ForMultipleMetric('V-Measure',algo_metric_dict.get('scg_vmeas'),algo_metric_dict.get('scclone_vmeas'),algo_metric_dict.get('bnpc_vmeas'),data_list,y_max,y_min,output,axes[1,1])
    lines = []
    labels = []

    print(" Figure axis ",fig.axes)
    #for ax in fig.axes[0]:
    Line, Label = fig.axes[0].get_legend_handles_labels()
    lines.extend(Line)
    labels.extend(Label)

    print(" Lines ",lines," Labels ",labels)
    fig.tight_layout()
    fig.legend(lines, labels, loc='lower left', ncol=3, fontsize=18)
    fig.subplots_adjust(left=0.1, bottom=0.1, right=1)
    fig.savefig('boxplots/'+output+'.png',bbox_inches="tight")

# def boxplot_ForMultipleMetric(metric, lgTree_data''', scclone_data, bnpc_data, rc_data, scite_data, data_list''', y_max, y_min, output, axis):
def boxplot_ForMultipleMetric(metric, lgTree_data, lace_data, data_list, y_max, y_min, output, axis):
    ticks = ['3','4','5']
    #plt.figure()
    #axis.set_title('Beta values', fontsize=15)

    #scg = axis.boxplot(scg_data, positions=np.array(range(len(scg_data)))*2.0-0.6, sym='', widths=0.2, whis=(0,100))
    #scclone = axis.boxplot(scclone_data, positions=np.array(range(len(scclone_data)))*2.0-0.4, sym='', widths=0.2, whis=(0,100))
    #bnpc = axis.boxplot(bnpc_data, positions=np.array(range(len(bnpc_data)))*2.0-0.2, sym='', widths=0.2, whis=(0,100))
    lgTree = axis.boxplot(lgTree_data, positions=np.array(range(len(lgTree_data)))*2.0+0.0, sym='', widths=0.2, whis=(0,100))
    lace = axis.boxplot(lace_data, positions=np.array(range(len(lace_data)))*2.0+0.2, sym='', widths=0.2, whis=(0,100))
    
    #siclonefit = axis.boxplot(siclonefit_data, positions=np.array(range(len(siclonefit_data)))*2.0+0.4, sym='', widths=0.2, whis=(0,100))
    #scg = plt.boxplot(scg_data, positions=[0.3], sym='', widths=0.2)
    #scclone = plt.boxplot(scclone_data, positions=[0.3], sym='', widths=0.2)
    #bnpc = plt.boxplot(bnpc_data, positions=[0.3], sym='', widths=0.2)

    set_box_color(lgTree, '#e41a1c') # colors are from http://colorbrewer2.org/
    set_box_color(lace, '#377eb8')
    #set_box_color(bnpc, '#4daf4a')
    #set_box_color(lgTree, '#feb24c')
    #set_box_color(scite, '#c51b8a')
    #set_box_color(siclonefit, '#636363')

    # draw temporary red and blue lines and use them to create a legend
    axis.plot([], c='#e41a1c', label='lgTree')
    axis.plot([], c='#377eb8', label='LACE')
    #axis.plot([], c='#4daf4a', label='BnpC')
    #axis.plot([], c='#feb24c', label='lgTree')
    #axis.plot([], c='#c51b8a', label='SCITE')
    #axis.plot([], c='#636363', label='SiCloneFit')
    #plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, fontsize=15)

    axis.set_xticks(range(0, len(ticks) * 2, 2), ticks, fontsize=18)
    axis.set_xlim(-2, len(ticks)*2)
    axis.tick_params(axis="y",labelsize='xx-large')
    #axis.set_yticks(fontsize=15)
    if metric == "V-Measure":
        axis.set_ylim(0.5, 1.01)
    elif metric == "Precision":
        plt.ylim(0, 110)
    elif metric == "Recall":
        axis.set_ylim(-0.10, 1.01)
    elif metric == "Time":
        plt.ylim(y_min, y_max+1000)
    #else:
    #    axis.set_ylim(50, 110)
    #if metric == 'Accuracy' or metric == 'Sensitivity' or metric == 'Specificity':
    #    axis.set_ylabel(metric+'(%)', fontsize=14)
    if metric == 'Precision' or metric == 'Recall':
        axis.set_ylabel(metric+' (%)', fontsize=18)
    elif metric == 'Time':
        axis.set_ylabel('Running time (s)', fontsize=18)
    else:
        axis.set_ylabel(metric, fontsize=18)
    #axis.tight_layout()

#def boxplot_ForMetric(metric, scg_data, scclone_data, bnpc_data, data_list, output):
#    ticks = ['0.1','0.3','0.4','0.2']
#    plt.figure()
#    plt.title('Beta values', fontsize=15)

#    scg = plt.boxplot(scg_data, positions=np.array(range(len(scg_data)))*2.0-0.4, sym='', widths=0.4)
#    scclone = plt.boxplot(scclone_data, positions=np.array(range(len(scclone_data)))*2.0+0.0, sym='', widths=0.4)
#    bnpc = plt.boxplot(bnpc_data, positions=np.array(range(len(bnpc_data)))*2.0+0.4, sym='', widths=0.4)
    #scg = plt.boxplot(scg_data, positions=[0.3], sym='', widths=0.2)
    #scclone = plt.boxplot(scclone_data, positions=[0.3], sym='', widths=0.2)
    #bnpc = plt.boxplot(bnpc_data, positions=[0.3], sym='', widths=0.2)

#    set_box_color(scg, '#e41a1c') # colors are from http://colorbrewer2.org/
#    set_box_color(scclone, '#377eb8')
#    set_box_color(bnpc, '#4daf4a')

    # draw temporary red and blue lines and use them to create a legend
#    plt.plot([], c='#e41a1c', label='SCG')
#    plt.plot([], c='#377eb8', label='SCClone')
#    plt.plot([], c='#4daf4a', label='BnpC')
    #plt.legend(bbox_to_anchor=(1.02, 1), loc='upper left', borderaxespad=0, fontsize=15)

#    plt.xticks(range(0, len(ticks) * 2, 2), ticks, fontsize=15)
#    plt.xlim(-2, len(ticks)*2)
#    plt.yticks(fontsize=15)
#    if metric == "V-Measure":
#        plt.ylim(0.2, 1.2)
    #elif metric == "Sensitivity":
    #    plt.ylim(20, 110)
#    else:
#        plt.ylim(50, 110)
#    if metric == 'Accuracy' or metric == 'Sensitivity' or metric == 'Specificity':
#        plt.ylabel(metric+'%', fontsize=15)
#    else:
#        plt.ylabel(metric, fontsize=15)
#    plt.tight_layout()
    #plt.savefig('metrics_boxplots/'+output+'_boxcompare_'+metric+'.png')

#def plot_boxplots(algo_metric_dict, data_list, output):
    #boxplot_ForMetric('Accuracy',algo_metric_dict.get('scg_acc'),algo_metric_dict.get('scclone_acc'),algo_metric_dict.get('bnpc_acc'),data_list,output)
#    boxplot_ForMetric('V-Measure',algo_metric_dict.get('scg_vmeas'),algo_metric_dict.get('scclone_vmeas'),algo_metric_dict.get('bnpc_vmeas'),data_list,output)
#    boxplot_ForMetric('Sensitivity',algo_metric_dict.get('scg_sens'),algo_metric_dict.get('scclone_sens'),algo_metric_dict.get('bnpc_sens'),data_list,output)
#    boxplot_ForMetric('Specificity',algo_metric_dict.get('scg_spec'),algo_metric_dict.get('scclone_spec'),algo_metric_dict.get('bnpc_spec'),data_list,output)

''' This method gets the metrics for each reptitions of each dataset in a dictionary for each algorithm passed as dictName. 
For example, For t8 the accuracy of one algo will be saved in the key 't8_acc': [rep1_val, .. , rep2_val] '''
def get_data(rootdir,dictName,data_list):
    print(" Rootdir ",rootdir)
    dictName = {}
    print(" Data list ",data_list)
    #sum_time = 0
    for subdir, dirs, files in os.walk(rootdir):
        #print(" subdir ",subdir)
        if algo == "LACE" and len(subdir.split("/")) > 3:
            algo_name = (subdir.split("/"))[1]
            sub_keyname = (subdir.split("/"))[3]
            #print(" Sub keyname ",sub_keyname)
        elif len(subdir.split("/")) > 3:
            algo_name = (subdir.split("/"))[0]
            sub_keyname = (subdir.split("/"))[2]
            #print(" Sub keyname ",sub_keyname)
        else:
            continue
        if sub_keyname not in data_list:
            #print(sub_keyname)
            continue
        for file in files:
            if "eval" in file:
                eval_path = os.path.join(subdir, file)
                #print(" Eval metrics path ",eval_path)
                precision = 0
                recall = 0
                rounded_vmeas_val = 0
                with open(eval_path,'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if "Total Precision" in line:
                            line_list = line.split(" ")
                            print(" Line list ",line_list)
                            precision = float(line_list[4].strip())
                        if "Total Recall" in line:
                            recall = float(line_list[4].strip())
                            #print(" Recall value ",recall)
                        if "V measure" in line:
                            vmeas_val = ((line.split(" "))[5]).strip()
                            rounded_vmeas_val = round(float(vmeas_val), 3)

                #mean_precision = prec_sum/tp_count
                #mean_recall = recall_sum/tp_count
                prec_key_name = sub_keyname+"_prec"
                recall_key_name = sub_keyname+"_recall"
                vmeas_key_name = sub_keyname+"_vmeas"

                if algo == "LACE" and sub_keyname == "t5":
                    if prec_key_name in dictName:
                        temp_list = dictName[prec_key_name]
                        temp_list.append(0)
                        dictName[prec_key_name] = temp_list
                    else:
                        prec_list = []
                        prec_list.append(0)
                        dictName[prec_key_name] = prec_list

                    if recall_key_name in dictName:
                        temp_list = dictName[recall_key_name]
                        temp_list.append(0)
                        dictName[recall_key_name] = temp_list
                    else:
                        recall_list = []
                        recall_list.append(0)
                        dictName[recall_key_name] = recall_list

                    if vmeas_key_name in dictName:
                        temp_list = dictName[vmeas_key_name]
                        temp_list.append(0)
                        dictName[vmeas_key_name] = temp_list
                    else:
                        vmeas_list = []
                        vmeas_list.append(0)
                        dictName[vmeas_key_name] = vmeas_list
                else:
                    if prec_key_name in dictName:
                        temp_list = dictName[prec_key_name]
                        temp_list.append(precision)
                        dictName[prec_key_name] = temp_list
                    else:
                        prec_list = []
                        prec_list.append(precision)
                        dictName[prec_key_name] = prec_list

                    if recall_key_name in dictName:
                        temp_list = dictName[recall_key_name]
                        temp_list.append(recall)
                        dictName[recall_key_name] = temp_list
                    else:
                        recall_list = []
                        recall_list.append(recall)
                        dictName[recall_key_name] = recall_list

                    if vmeas_key_name in dictName:
                        temp_list = dictName[vmeas_key_name]
                        temp_list.append(rounded_vmeas_val)
                        dictName[vmeas_key_name] = temp_list
                    else:
                        vmeas_list = []
                        vmeas_list.append(rounded_vmeas_val)
                        dictName[vmeas_key_name] = vmeas_list
            
            if algo == "lgTree" and "total_time" in file:
                time_path = os.path.join(subdir, file)
                #print(" Running time path ",time_path)

                with open(time_path,'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        #print(" Line ",line)
                        time_obj = int(line)
                
                key_name = sub_keyname+"_time"
                #key_name = sub_keyname
                if key_name in dictName:
                    temp_list = dictName[key_name]
                    temp_list.append(time_obj)
                    dictName[key_name] = temp_list
                else:
                    time_val_list = []
                    time_val_list.append(time_obj)
                    dictName[key_name] = time_val_list

            if algo == "LACE" and "time" in file:
                time_path = os.path.join(subdir, file)
                print(" LACE Running time path ",time_path)

                with open(time_path,'r') as f:
                    lines = f.readlines()
                    for line in lines:
                        if "Job Wall-clock time" in line:
                            time_line = (line.split(" "))[3]
                            hour_time = time_line.split(":")[0]
                            if '-' in hour_time:
                                time_in_secs = -1
                            else:
                                hour_time = int((time_line.split(":"))[0])*60*60
                                min_time = int((time_line.split(":"))[1])*60
                                sec_time = int((time_line.split(":"))[2])
                                time_in_secs = hour_time+min_time+sec_time
                            if time_in_secs >= 86400:
                                time_in_secs = -1
                    key_name = sub_keyname+"_time"
                    if key_name in dictName:
                        temp_list = dictName[key_name]
                        temp_list.append(time_in_secs)
                        dictName[key_name] = temp_list
                    else:
                        time_val_list = []
                        time_val_list.append(time_in_secs)
                        dictName[key_name] = time_val_list

    print(" Algo Dict ",dictName)
    return dictName

''' Arrange the data in a list for boxplots. '''
def create_lists(algoDict,data_list):
    #acc_dict = {}
    vmeas_dict = {}
    prec_dict = {}
    recall_dict = {}
    time_dict = {}
    for k,v in algoDict.items():
        #if 'acc' in k:
        #    acc_dict[k] = v
        if 'prec' in k:
            prec_dict[k] = v
        elif 'recall' in k:
            recall_dict[k] = v
        elif 'time' in k:
            time_dict[k] = v
        elif 'vmeas' in k:
            vmeas_dict[k] = v

    print(" Time dict ",time_dict)
    #acc_list = []
    vmeas_list = []
    prec_list = []
    recall_list = []
    time_list = []
    max_time = 0
    min_time = 0
    for data in data_list:
        #acc_v = [val for key, val in acc_dict.items() if data in key]
        #acc_list.append(acc_v[0])
        vmeas_v = [val for key, val in vmeas_dict.items() if data == (key.split('_'))[0]]
        vmeas_list.append(vmeas_v[0])
        prec_v = [val for key, val in prec_dict.items() if data == (key.split('_'))[0]]
        prec_list.append(prec_v[0])
        recall_v = [val for key, val in recall_dict.items() if data == (key.split('_'))[0]]
        recall_list.append(recall_v[0])
        v = [val for key, val in time_dict.items() if data == (key.split('_'))[0]]
        print(v)
        #max_time = 21000
        #min_time = 0
        if max_time < max(v[0]):
            max_time = max(v[0])
        if min_time > min(v[0]):
            min_time = min(v[0])
        time_list.append(v[0])
    return prec_list, recall_list, vmeas_list, time_list, max_time, min_time

parser = argparse.ArgumentParser()
parser.add_argument("-algo", "--algo",dest ="algo", nargs="*", help="Algorithm name to use")
parser.add_argument("-datasets", "--datasets",dest ="datasets", nargs="*", help="Datasets for which to get the metrics")
parser.add_argument("-output", "--output",dest ="output", help="Annotate the graphs for each condition")
#parser.add_argument("-title", "--title",dest ="title", help="Algorithm name to use")
args = parser.parse_args()

data_list = args.datasets
algo_list = args.algo
#print(args.algo)
#print(args.datasets)

# Loop over algo list and save each lists into data_lists for the box_plot. Then think about how to get the box plots for each metric.
algo_metric_dict = {}
max_min_list = []
for algo in algo_list:
    print(algo)
    if algo == "LACE":
        algo_dict = get_data('../LACE_exp/sim_output',algo,data_list)
    #elif algo == "siclonefit":
    #    algo_dict = get_data(algo+'/less24hrs/sim_output',algo,data_list)
    #elif algo == "SCITE":
    #    algo_dict = get_data(algo+'/sim_output',algo,data_list)
    else:
        algo_dict = get_data('./sim_output',algo,data_list) # One more variable for the algorithm 
    # For default bnpc add a condition here to access the default sim_output folder
    prec_list, recall_list, vmeas_list, time_list, max_time, min_time = create_lists(algo_dict,data_list)
    #algo_metric_dict[algo+"_acc"] = acc_list
    algo_metric_dict[algo+"_vmeas"] = vmeas_list
    algo_metric_dict[algo+"_prec"] = prec_list
    algo_metric_dict[algo+"_recall"] = recall_list
    algo_metric_dict[algo+"_time"] = time_list
    max_min_list.append(max_time)
    max_min_list.append(min_time)
    #print(" Accuracy ",acc_list)
    #print(" Sensitivity ",sens_list)
    #print(" Specificity ",spec_list)
    #print(" V measure ",vmeas_list)

y_max = max(max_min_list)
y_min = min(max_min_list)
print(" Y max ",y_max," Y min ",y_min)
print("====================================")
print(pd.DataFrame.from_dict(algo_metric_dict))
#algo_metric_df = pd.DataFrame.from_dict(algo_metric_dict)
#algo_metric_df.to_csv('doublets.csv')
print(algo_metric_dict)
#plot_boxplots(algo_metric_dict, data_list, args.output)
plot_multipleBp(algo_metric_dict, data_list, y_max, y_min, args.output)
