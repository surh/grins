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
import pandas as pd

sns.set(style="whitegrid",font_scale=1.5)
font = {'family': 'sans-serif',
        'color':  'black',
        'weight': 'normal',
        'size': 20, 'rotation':0,
        'verticalalignment': 'bottom',
        'horizontalalignment': 'left'
        }
# Minimal region of non-interrupted skew
region_length=250

# Threshold for AUC-to-length ratio:
threshold_AUC=0.1

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
def x_crossing(GCskew_array,TAskew_array,step,region_length):
    # Defining the threshold (e.g. 500/step means that we detect intervals of 500bp or more)
    interval_threshold=np.round(region_length/step,0)
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
def save_fig_overlap(record_name,GCskew_array,TAskew_array,overlapping_pos_intervals,overlapping_neg_intervals,parameter_set):
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
    plt.savefig("./Figs_overlapping_regions/"+record_name+"_found_skew_intervals_window"+str(parameter_set[0])+"_step"+str(parameter_set[1])+".png")
    plt.close()

# Plotting figure of predicted GRINS regions, according to AUCs
def save_fig_AUC(record_name,GCskew_array,TAskew_array,predicted_GRINS,parameter_set):
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
    plt.savefig("./Figs_detected_GRINS/"+record_name+"_found_GRINS_window"+str(parameter_set[0])+"_step"+str(parameter_set[1])+".png")
    plt.close()


# Adding values to bar graphs
def autolabel(series):
    for value in series:
        height = np.round(value.get_height(),2)
        ax.text(value.get_x() + value.get_width()/2., 1.05*height,
                str(height),
                ha='center', va='bottom')


if __name__ == "__main__":
    # Making a list of annotated GenBank files
    mypath = "input/"
    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    # print files
    # files_test=files[0:2]

    # Setting parameters for calculating skews
    parameters=[(150,30),(50,10),(30,30),(300,50)]


    # Vectors to keep results
    False_pos=[0 for a in range(0,len(parameters))]
    False_neg=[0 for a in range(0,len(parameters))]
    True_pos=[0 for a in range(0,len(parameters))]
    True_neg=[0 for a in range(0,len(parameters))]


    # Calculating skewsand plotting them; calculating AUCs and plotting them
    print files
    for item in files:
        print "Working on file",item
        record = SeqIO.read(mypath+item, "gb")
        input_DNA=record.seq
        record_name=record.name
        feature_list=record.features
        n=-1

        for parameter_set in parameters:
            n+=1
            window=parameter_set[0]
            step=parameter_set[1]
            GCskew_array=[]
            TAskew_array=[]
            ATskew_array=[]

            #Calculating skews
            for i in range(0,len(input_DNA),step):
                fragment=input_DNA[i:i+window]
                GCskew_array.append(GC_skew(fragment))
                TAskew_array.append(TA_skew(fragment))

            # Finding regions in which to calculate AUCs:
            # regions of >=500bp where both GC and TA skews are of the same sign (either both positive, or both negative)
            overlapping_pos_intervals,overlapping_neg_intervals=x_crossing(GCskew_array,TAskew_array,step,region_length)
            # save_fig_overlap(record_name,GCskew_array,TAskew_array,overlapping_pos_intervals,overlapping_neg_intervals,parameter_set)

            # Finding regions where AUC to interval length is larger than a threshold
            predicted_GRINS=[]
            # print parameter_set
            for i in range(0,len(overlapping_pos_intervals)):
                interval_length=float(overlapping_pos_intervals[i][1]-overlapping_pos_intervals[i][0])
                x=range(overlapping_pos_intervals[i][0],overlapping_pos_intervals[i][1])
                y_GC=GCskew_array[overlapping_pos_intervals[i][0]:overlapping_pos_intervals[i][1]]
                y_TA=TAskew_array[overlapping_pos_intervals[i][0]:overlapping_pos_intervals[i][1]]
                AUC_GCskew=metrics.auc(x,y_GC)
                AUC_TAskew=metrics.auc(x,y_TA)
                if AUC_GCskew/interval_length>threshold_AUC and AUC_TAskew/interval_length>threshold_AUC:
                    predicted_GRINS.append(overlapping_pos_intervals[i])
            for i in range(0,len(overlapping_neg_intervals)):
                interval_length=float(overlapping_neg_intervals[i][1]-overlapping_neg_intervals[i][0])
                x=range(overlapping_neg_intervals[i][0],overlapping_neg_intervals[i][1])
                y_GC=GCskew_array[overlapping_neg_intervals[i][0]:overlapping_neg_intervals[i][1]]
                y_TA=TAskew_array[overlapping_neg_intervals[i][0]:overlapping_neg_intervals[i][1]]
                AUC_GCskew=metrics.auc(x,y_GC)
                AUC_TAskew=metrics.auc(x,y_TA)
                if AUC_GCskew/interval_length>threshold_AUC and AUC_TAskew/interval_length>threshold_AUC:
                    predicted_GRINS.append(overlapping_neg_intervals[i])

            # Showing AUCs
            # save_fig_AUC(record_name,GCskew_array,TAskew_array,predicted_GRINS,parameter_set)

            #  Generating a list of regions with annotated GRINS
            annotated_GRINS=[]
            for feature in feature_list:
                if feature.type=="GRINS":
                    annotated_GRINS.append((np.round(feature.location.start/30,0),np.round(feature.location.end/30,0)))

            # Comparing GRINS annotations with detected GRINSs, and determining FP, FN, TP, TN
            FP=0
            FN=0
            TP=0
            TN=0
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
            False_pos[n]+=FP
            False_neg[n]+=FN
            True_pos[n]+=TP
            True_neg[n]+=TN
    # print False_pos
    # print False_neg
    # print True_pos
    # print True_neg


    # Making figures of sensitivity and selectivity
    sensitivities=[0 for a in range(0,len(parameters))]
    specificities=[0 for a in range(0,len(parameters))]
    precision=[0 for a in range(0,len(parameters))]
    for i in range(0,len(True_pos)):
        sensitivities[i]=(float(True_pos[i])/float(True_pos[i]+False_neg[i]))
        specificities[i]=(float(True_neg[i])/float(True_neg[i]+False_pos[i]))
        precision[i]=(float(True_pos[i])/float(True_pos[i]+False_pos[i]))
    recall=sensitivities

    with open("Detection_results.txt",'w') as results_file:
        results_file.write("Detection with region threshold of "+str(region_length)+"bp and AUC-to-length threshold of "+str(threshold_AUC)+"\n")
        results_file.write("for parameter sets \n")
        for k in range(0,len(parameters)):
            results_file.write("\n")
            results_file.write("window "+str(parameters[k][0])+" , step "+str(parameters[k][1])+":\n")
            results_file.write("sensitivity "+str(np.round(sensitivities[k],3))+"\n")
            results_file.write("specificity "+str(np.round(specificities[k],3))+"\n")
            results_file.write("precision "+str(np.round(precision[k],3))+"\n")


    x=[1,2,3,4]
    plt.figure(1)
    ax=plt.subplot(111)
    ax.set_ylim(0,1.2)
    sensit=ax.bar([a-0.15 for a in x],sensitivities,width=0.2,color="red", align="center")
    autolabel(sensit)
    specif=ax.bar([a+0.15 for a in x],specificities,width=0.2,color="blue", align="center")
    autolabel(specif)
    ax.set_xticks(x)
    ax.set_xticklabels(["w150\ns30","w50\ns10","w30\ns30","w300\ns50"])
    ax.legend(["Sensitivity","Specificity"])
    ax.set_title("Region threshold = "+str(region_length)+"bp\nAUC-to-length threshold = "+str(threshold_AUC))
    plt.savefig("Detection_results_sens_spec.pdf")
    plt.close()

    plt.figure(2)
    ax=plt.subplot(111)
    ax.set_ylim(0,1.2)
    prec=ax.bar([a-0.15 for a in x],precision,width=0.2,color="purple", align="center")
    autolabel(prec)
    rec=ax.bar([a+0.15 for a in x],recall,width=0.2,color="red", align="center")
    autolabel(rec)
    ax.set_xticks(x)
    ax.set_xticklabels(["w150\ns30","w50\ns10","w30\ns30","w300\ns50"])
    ax.legend(["Precision","Recall"])
    ax.set_title("Region threshold = "+str(region_length)+"bp\nAUC-to-length threshold = "+str(threshold_AUC))
    plt.savefig("Detection_results_prec_rec.pdf")
    plt.close()
