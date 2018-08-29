#!/usr/bin/env python
import os
from os import listdir
from os.path import isfile, join
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as clr
import seaborn as sns
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature,FeatureLocation
from sklearn import metrics
import argparse

sns.set(style="whitegrid",font_scale=1.5)
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 20, 'rotation':0,
        'verticalalignment': 'bottom',
        'horizontalalignment': 'left'
        }


# Calculating skews
def GC_skew(string):
    G_number=string.count("G")
    C_number=string.count("C")
    if G_number+C_number!=0:
        skew=float(G_number-C_number)/float(G_number+C_number)
    else:
        skew=0
    return skew

def TA_skew(string):
    A_number=string.count("A")
    T_number=string.count("T")
    if A_number+T_number!=0:
        skew=float(T_number-A_number)/float(A_number+T_number)
    else:
        skew=0
    return skew


# Finding places where each skew crosses the x axis, and finding regions of overlap
def x_crossing(GCskew_array,TAskew_array,step,length_threshold):
    # Defining the threshold (e.g. 500/step means that we detect intervals of 500bp or more)
    interval_threshold=np.round(length_threshold/step,0)
    # Finding positions where GC skew values cross the x axis
    GCskew_x_crossing_pos=[]
    GCskew_x_crossing_neg=[]
    for i in range(0,len(GCskew_array)-1):
        if (GCskew_array[i]>=0 and GCskew_array[i+1]<0):
            GCskew_x_crossing_neg.append(i+1)
        elif (GCskew_array[i]<0 and GCskew_array[i+1]>=0):
            GCskew_x_crossing_pos.append(i+1)
    # Finding intervals of 500bp or longer, where the GC skew sign does not change
    significant_GCskew_intervals_pos=[]
    significant_GCskew_intervals_neg=[]
    # case when first crossing is "pos"
    if GCskew_x_crossing_pos[0]<GCskew_x_crossing_neg[0]:
        for j in range(0,len(GCskew_x_crossing_neg)):
            if GCskew_x_crossing_neg[j]-GCskew_x_crossing_pos[j]>=interval_threshold:
                significant_GCskew_intervals_pos.append((GCskew_x_crossing_pos[j],GCskew_x_crossing_neg[j]))
        for j in range(1,len(GCskew_x_crossing_pos)):
            if GCskew_x_crossing_pos[j]-GCskew_x_crossing_neg[j-1]>=interval_threshold:
                significant_GCskew_intervals_neg.append((GCskew_x_crossing_neg[j-1],GCskew_x_crossing_pos[j]))
    # case when first crossing is "neg"
    elif GCskew_x_crossing_neg[0]<GCskew_x_crossing_pos[0]:
        for j in range(0,len(GCskew_x_crossing_pos)):
            if GCskew_x_crossing_pos[j]-GCskew_x_crossing_neg[j]>=interval_threshold:
                significant_GCskew_intervals_neg.append((GCskew_x_crossing_neg[j],GCskew_x_crossing_pos[j]))
        for j in range(1,len(GCskew_x_crossing_neg)):
            if GCskew_x_crossing_neg[j]-GCskew_x_crossing_pos[j-1]>=interval_threshold:
                significant_GCskew_intervals_pos.append((GCskew_x_crossing_pos[j-1],GCskew_x_crossing_neg[j]))
    # Finding positions where TA skew values cross the x axis
    TAskew_x_crossing_pos=[]
    TAskew_x_crossing_neg=[]
    for i in range(0,len(TAskew_array)-1):
        if (TAskew_array[i]>=0 and TAskew_array[i+1]<0):
            TAskew_x_crossing_neg.append(i+1)
        elif (TAskew_array[i]<0 and TAskew_array[i+1]>=0):
            TAskew_x_crossing_pos.append(i+1)
    # Finding intervals of 500bp or longer, where the GC skew sign does not change
    significant_TAskew_intervals_pos=[]
    significant_TAskew_intervals_neg=[]
    # case when first crossing is "pos"
    if TAskew_x_crossing_pos[0]<TAskew_x_crossing_neg[0]:
        for j in range(0,len(TAskew_x_crossing_neg)):
            if TAskew_x_crossing_neg[j]-TAskew_x_crossing_pos[j]>=interval_threshold:
                significant_TAskew_intervals_pos.append((TAskew_x_crossing_pos[j],TAskew_x_crossing_neg[j]))
        for j in range(1,len(TAskew_x_crossing_pos)):
            if TAskew_x_crossing_pos[j]-TAskew_x_crossing_neg[j-1]>=interval_threshold:
                significant_TAskew_intervals_neg.append((TAskew_x_crossing_neg[j-1],TAskew_x_crossing_pos[j]))
    # case when first crossing is "neg"
    elif TAskew_x_crossing_neg[0]<TAskew_x_crossing_pos[0]:
        for j in range(0,len(TAskew_x_crossing_pos)):
            if TAskew_x_crossing_pos[j]-TAskew_x_crossing_neg[j]>=interval_threshold:
                significant_TAskew_intervals_neg.append((TAskew_x_crossing_neg[j],TAskew_x_crossing_pos[j]))
        for j in range(1,len(TAskew_x_crossing_neg)):
            if TAskew_x_crossing_neg[j]-TAskew_x_crossing_pos[j-1]>=interval_threshold:
                significant_TAskew_intervals_pos.append((TAskew_x_crossing_pos[j-1],TAskew_x_crossing_neg[j]))
    # Finding intervals where GC and TA skews are both of positive sign, and overlap over 500bp or longer
    overlapping_pos_intervals=[]
    for k1 in range(0,len(significant_GCskew_intervals_pos)):
        for k2 in range(0,len(significant_TAskew_intervals_pos)):
            # case where the start of TA region is within the GC region
            if (significant_GCskew_intervals_pos[k1][0]<=significant_TAskew_intervals_pos[k2][0]) and (significant_GCskew_intervals_pos[k1][1]>significant_TAskew_intervals_pos[k2][0]):
                interval_start=significant_TAskew_intervals_pos[k2][0]
                interval_end=min(significant_GCskew_intervals_pos[k1][1],significant_TAskew_intervals_pos[k2][1])
                if interval_end-interval_start>=interval_threshold:
                    overlapping_pos_intervals.append((interval_start,interval_end))
            #case where the start of GC region is within the TA region
            elif (significant_TAskew_intervals_pos[k2][0]<=significant_GCskew_intervals_pos[k1][0]) and (significant_TAskew_intervals_pos[k2][1]>significant_GCskew_intervals_pos[k1][0]):
                interval_start=significant_GCskew_intervals_pos[k1][0]
                interval_end=min(significant_TAskew_intervals_pos[k2][1],significant_GCskew_intervals_pos[k1][1])
                if interval_end-interval_start>=interval_threshold:
                    overlapping_pos_intervals.append((interval_start,interval_end))
    # Finding intervals where GC and TA skews are both of negative sign, and overlap over 500bp or longer
    overlapping_neg_intervals=[]
    for k1 in range(0,len(significant_GCskew_intervals_neg)):
        for k2 in range(0,len(significant_TAskew_intervals_neg)):
            # case where the start of TA region is withon the GC region
            if (significant_GCskew_intervals_neg[k1][0]<=significant_TAskew_intervals_neg[k2][0]) and (significant_GCskew_intervals_neg[k1][1]>significant_TAskew_intervals_neg[k2][0]):
                interval_start=significant_TAskew_intervals_neg[k2][0]
                interval_end=min(significant_GCskew_intervals_neg[k1][1],significant_TAskew_intervals_neg[k2][1])
                if interval_end-interval_start>=interval_threshold:
                    overlapping_neg_intervals.append((interval_start,interval_end))
            #case where the start of GC region is withon the TA region
            elif (significant_TAskew_intervals_neg[k2][0]<=significant_GCskew_intervals_neg[k1][0]) and (significant_TAskew_intervals_neg[k2][1]>significant_GCskew_intervals_neg[k1][0]):
                interval_start=significant_GCskew_intervals_neg[k1][0]
                interval_end=min(significant_TAskew_intervals_neg[k2][1],significant_GCskew_intervals_neg[k1][1])
                if interval_end-interval_start>=interval_threshold:
                    overlapping_neg_intervals.append((interval_start,interval_end))
    return overlapping_pos_intervals,overlapping_neg_intervals


# Plotting figure of overlapping regions
def save_fig_overlap(folder_name,record_name,GCskew_array,TAskew_array,overlapping_pos_intervals,overlapping_neg_intervals):
    a1=range(0,len(GCskew_array))
    plt.figure(1)
    fig,ax1=plt.subplots(figsize=(10, 4))
    ax1.set_ylim(-1,1)
    ax1.set_xlim(0,len(GCskew_array))
    ax1.grid(False)
    plt.plot(a1,GCskew_array,color="blue",linewidth=0.5)
    plt.plot(a1,TAskew_array,color="green",linewidth=0.5)
    plt.plot([0,len(GCskew_array)],[0,0],color="black",linewidth=0.5)
    for i in range(0,len(overlapping_pos_intervals)):
        rect = patches.Rectangle((overlapping_pos_intervals[i][0],1),overlapping_pos_intervals[i][1]-overlapping_pos_intervals[i][0],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
        ax1.add_patch(rect)
    for i in range(0,len(overlapping_neg_intervals)):
        rect = patches.Rectangle((overlapping_neg_intervals[i][0],1),overlapping_neg_intervals[i][1]-overlapping_neg_intervals[i][0],-2,linewidth=1,edgecolor='pink',facecolor='pink', alpha=0.3)
        ax1.add_patch(rect)
    plt.legend(["GC skew","TA skew"],loc=3)
    plt.savefig(folder_name+"/Figures/"+record_name+"_found_intervals.png")
    plt.close()

# Plotting figure of predicted GRINS regions, according to AUCs
def save_fig_AUC(folder_name,record_name,GCskew_array,TAskew_array,predicted_GRINS):
    a1=range(0,len(GCskew_array))
    plt.figure(1)
    fig,ax1=plt.subplots(figsize=(10, 4))
    ax1.set_ylim(-1,1)
    ax1.set_xlim(0,len(GCskew_array))
    ax1.grid(False)
    plt.plot(a1,GCskew_array,color="blue",linewidth=0.5)
    plt.plot(a1,TAskew_array,color="green",linewidth=0.5)
    plt.plot([0,len(GCskew_array)],[0,0],color="black",linewidth=0.5)
    for i in range(0,len(predicted_GRINS)):
        rect = patches.Rectangle((predicted_GRINS[i][0],1),predicted_GRINS[i][1]-predicted_GRINS[i][0],-2,linewidth=1,edgecolor='gray',facecolor='gray', alpha=0.3)
        ax1.add_patch(rect)
    plt.legend(["GC skew","TA skew"],loc=3)
    plt.savefig(folder_name+"/Figures/"+record_name+"_found_GRINSs.png")
    plt.close()


# Adding values to bar graphs
def autolabel(series):
    for value in series:
        height = np.round(value.get_height(),2)
        ax.text(value.get_x() + value.get_width()/2., 1.05*height,
                str(height),
                ha='center', va='bottom')

def AUC_cutoff(overlapping_interval,AUC_threshold):
    global predicted_GRINS
    interval_length=float(overlapping_interval[1]-overlapping_interval[0])
    x=range(overlapping_interval[0],overlapping_interval[1])
    y_GC=GCskew_array[overlapping_interval[0]:overlapping_interval[1]]
    y_TA=TAskew_array[overlapping_interval[0]:overlapping_interval[1]]
    AUC_GCskew=abs(metrics.auc(x,y_GC))
    AUC_TAskew=abs(metrics.auc(x,y_TA))
    if AUC_GCskew/interval_length>=AUC_threshold and AUC_TAskew/interval_length>=AUC_threshold:
        predicted_GRINS.append(overlapping_interval)


if __name__ == "__main__":
    # Making a list of annotated GenBank files
    mypath = "input/"
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    # print files
    # files_test=files[0:2]

    # Setting parameters 
    parser=argparse.ArgumentParser(description="User-specified parameters for GRINS detection")
    # Size of the sliding window
    parser.add_argument('-w',"--window",type=int, help="window size for calculating skews", default=150)
    # Size of the step
    parser.add_argument('-s',"--step",type=int, help="step size for calculating skews", default=30)
    # Threshold for minimal region of non-interrupted skew
    parser.add_argument('-lent',"--length_threshold",type=float, help="threshold for fragment length", default=250)
    # Threshold for AUC-to-length ratio:
    parser.add_argument('-AUCt ',"--AUC_threshold",type=float, help="threshold for AUC-to-length ratio", default=0.1)
    # Folder name:
    parser.add_argument('-name ',"--folder_name",type=str, default="Results", help="name of the folder to store files")

    args=parser.parse_args()
    
    if args.folder_name=="Results":
        folder_name="./Window{}_step{}_length{}_AUC{}/".format(args.window,args.step,args.length_threshold,args.AUC_threshold)
    else:
        folder_name="./{}/".format(args.folder_name)
    os.mkdir(folder_name)
    os.mkdir(folder_name+"/Figures/")

    # Variables and vectors to keep results
    FP=0
    FN=0
    TP=0
    TN=0



    with open(folder_name+"Detection_results.txt",'w') as results_file:
        results_file.write("Detection with sliding window of "+str(args.window)+" bp, step of "+str(args.step)+", region threshold of "+str(args.length_threshold)+" bp and AUC-to-length threshold of "+str(args.AUC_threshold)+"\n")
        
        # Calculating skews and plotting them; calculating AUCs and plotting them
        print(files)
        for item in files:
            print("Working on file",item)
            record = SeqIO.read(mypath+item, "gb")
            input_DNA=record.seq
            record_name=record.name
            feature_list=record.features
            results_file.write(record_name+"\n")
            results_file.write("Detected GRINSs:\n")

            #Calculating skews
            GCskew_array=[]
            TAskew_array=[]
            ATskew_array=[]
            for i in range(0,len(input_DNA),args.step):
                fragment=input_DNA[i:i+args.window]
                GCskew_array.append(GC_skew(fragment))
                TAskew_array.append(TA_skew(fragment))

            # Finding regions in which to calculate AUCs:
            # regions of >=500bp where both GC and TA skews are of the same sign (either both positive, or both negative)
            overlapping_pos_intervals,overlapping_neg_intervals=x_crossing(GCskew_array,TAskew_array,args.step,args.length_threshold)
            save_fig_overlap(folder_name,record_name,GCskew_array,TAskew_array,overlapping_pos_intervals,overlapping_neg_intervals)

            # Finding regions where AUC to interval length is larger than a threshold
            predicted_GRINS=[]
            # print parameter_set
            for i in range(0,len(overlapping_pos_intervals)):
                AUC_cutoff(overlapping_pos_intervals[i],args.AUC_threshold)
            for i in range(0,len(overlapping_neg_intervals)):
                AUC_cutoff(overlapping_neg_intervals[i],args.AUC_threshold)

            # Showing AUCs
            save_fig_AUC(folder_name,record_name,GCskew_array,TAskew_array,predicted_GRINS)
            
            # Recording annotated GRINSs
            for item in predicted_GRINS:
                results_file.write(str(item[0])+"-"+str(item[1])+"\n")


            #  Generating a list of regions with annotated GRINS
            annotated_GRINS=[]
            
            results_file.write("Annotated GRINSs:\n")
            for feature in feature_list:
                if feature.type=="GRINS":
                    annotated_GRINS.append((np.round(feature.location.start/30,0),np.round(feature.location.end/30,0)))
                    results_file.write(str(np.round(feature.location.start/30,0))+"-"+str(np.round(feature.location.end/30,0))+"\n")
            results_file.write("\n")

            # Comparing GRINS annotations with detected GRINSs, and determining FP, FN, TP, TN
            for i in range(0,len(input_DNA)):
                true_GRINS=0
                annot_GRINS=0
                for k1 in range(0,len(annotated_GRINS)):
                    if i>=annotated_GRINS[k1][0] and i<annotated_GRINS[k1][1]:
                        true_GRINS=1
                        # print "pos in true GRINS",i
                        break
                for k2 in range(0,len(predicted_GRINS)):
                    if i>=predicted_GRINS[k2][0] and i<predicted_GRINS[k2][1]:
                        annot_GRINS=1
                        # print "pos in annot GRINS",i
                        break
                if true_GRINS==1 and annot_GRINS==1:
                    TP+=1
                elif true_GRINS==1 and annot_GRINS==0:
                    FN+=1
                elif true_GRINS==0 and annot_GRINS==1:
                    FP+=1
                else:
                    TN+=1


        # Saving results
        sensitivity=0
        specificity=0
        precision=0
        sensitivity=(float(TP)/float(TP+FN))
        specificity=(float(TN)/float(TN+FP))
        precision=(float(TP)/float(TP+FP))
        recall=sensitivity

    with open(folder_name+"Detection_result_summary.txt",'w') as results_file2:
        results_file2.write("Detection with sliding window of "+str(args.window)+" bp, step of "+str(args.step)+", region threshold of "+str(args.length_threshold)+" bp and AUC-to-length threshold of "+str(args.AUC_threshold)+"\n")
        results_file2.write("sensitivity "+str(np.round(sensitivity,3))+"\n")
        results_file2.write("specificity "+str(np.round(specificity,3))+"\n")
        results_file2.write("precision "+str(np.round(precision,3))+"\n")
    
