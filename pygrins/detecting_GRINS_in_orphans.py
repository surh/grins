import os
from os import listdir
from os.path import isfile, join
import numpy as np
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature,FeatureLocation
import argparse
from sklearn import metrics


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
def x_crossing(GCskew_array,TAskew_array):

    # Finding positions where GC skew values cross the x axis
    GCskew_x_crossing_pos=[]
    GCskew_x_crossing_neg=[]

    for i in range(0,len(GCskew_array)-1):
        if (GCskew_array[i]>=0 and GCskew_array[i+1]<0):
            GCskew_x_crossing_neg.append(i+1)
        elif (GCskew_array[i]<0 and GCskew_array[i+1]>=0):
            GCskew_x_crossing_pos.append(i+1)

    if len(GCskew_x_crossing_pos)<1 or len(GCskew_x_crossing_neg)<1:
        return [],[]


    # Finding intervals, where the GC skew sign does not change
    significant_GCskew_intervals_pos=[]
    significant_GCskew_intervals_neg=[]

    # case when first crossing is "pos"
    if GCskew_x_crossing_pos[0]<GCskew_x_crossing_neg[0]:
        for j in range(0,len(GCskew_x_crossing_neg)):
            # if GCskew_x_crossing_neg[j]-GCskew_x_crossing_pos[j]>=interval_threshold:
            significant_GCskew_intervals_pos.append((GCskew_x_crossing_pos[j],GCskew_x_crossing_neg[j]))
        for j in range(1,len(GCskew_x_crossing_pos)):
            # if GCskew_x_crossing_pos[j]-GCskew_x_crossing_neg[j-1]>=interval_threshold:
            significant_GCskew_intervals_neg.append((GCskew_x_crossing_neg[j-1],GCskew_x_crossing_pos[j]))
    
    # case when first crossing is "neg"
    elif GCskew_x_crossing_neg[0]<GCskew_x_crossing_pos[0]:
        for j in range(0,len(GCskew_x_crossing_pos)):
            # if GCskew_x_crossing_pos[j]-GCskew_x_crossing_neg[j]>=interval_threshold:
            significant_GCskew_intervals_neg.append((GCskew_x_crossing_neg[j],GCskew_x_crossing_pos[j]))
        for j in range(1,len(GCskew_x_crossing_neg)):
            # if GCskew_x_crossing_neg[j]-GCskew_x_crossing_pos[j-1]>=interval_threshold:
            significant_GCskew_intervals_pos.append((GCskew_x_crossing_pos[j-1],GCskew_x_crossing_neg[j]))

    # Finding positions where TA skew values cross the x axis
    TAskew_x_crossing_pos=[]
    TAskew_x_crossing_neg=[]
    for i in range(0,len(TAskew_array)-1):
        if (TAskew_array[i]>=0 and TAskew_array[i+1]<0):
            TAskew_x_crossing_neg.append(i+1)
        elif (TAskew_array[i]<0 and TAskew_array[i+1]>=0):
            TAskew_x_crossing_pos.append(i+1)

    # Finding intervals, where the TA skew sign does not change
    significant_TAskew_intervals_pos=[]
    significant_TAskew_intervals_neg=[]

    # case when first crossing is "pos"
    if TAskew_x_crossing_pos[0]<TAskew_x_crossing_neg[0]:
        for j in range(0,len(TAskew_x_crossing_neg)):
            # if TAskew_x_crossing_neg[j]-TAskew_x_crossing_pos[j]>=interval_threshold:
            significant_TAskew_intervals_pos.append((TAskew_x_crossing_pos[j],TAskew_x_crossing_neg[j]))
        for j in range(1,len(TAskew_x_crossing_pos)):
            # if TAskew_x_crossing_pos[j]-TAskew_x_crossing_neg[j-1]>=interval_threshold:
            significant_TAskew_intervals_neg.append((TAskew_x_crossing_neg[j-1],TAskew_x_crossing_pos[j]))
    
    # case when first crossing is "neg"
    elif TAskew_x_crossing_neg[0]<TAskew_x_crossing_pos[0]:
        for j in range(0,len(TAskew_x_crossing_pos)):
            # if TAskew_x_crossing_pos[j]-TAskew_x_crossing_neg[j]>=interval_threshold:
            significant_TAskew_intervals_neg.append((TAskew_x_crossing_neg[j],TAskew_x_crossing_pos[j]))
        for j in range(1,len(TAskew_x_crossing_neg)):
            # if TAskew_x_crossing_neg[j]-TAskew_x_crossing_pos[j-1]>=interval_threshold:
            significant_TAskew_intervals_pos.append((TAskew_x_crossing_pos[j-1],TAskew_x_crossing_neg[j]))
    
    # Finding intervals where GC and TA skews are both of positive sign, and overlap
    overlapping_pos_intervals=[]
    for k1 in range(0,len(significant_GCskew_intervals_pos)):
        for k2 in range(0,len(significant_TAskew_intervals_pos)):
            # case where the start of TA region is within the GC region
            if (significant_GCskew_intervals_pos[k1][0]<=significant_TAskew_intervals_pos[k2][0]) and (significant_GCskew_intervals_pos[k1][1]>significant_TAskew_intervals_pos[k2][0]):
                interval_start=significant_TAskew_intervals_pos[k2][0]
                interval_end=min(significant_GCskew_intervals_pos[k1][1],significant_TAskew_intervals_pos[k2][1])
                # if interval_end-interval_start>=interval_threshold:
                overlapping_pos_intervals.append((interval_start,interval_end))
            #case where the start of GC region is within the TA region
            elif (significant_TAskew_intervals_pos[k2][0]<=significant_GCskew_intervals_pos[k1][0]) and (significant_TAskew_intervals_pos[k2][1]>significant_GCskew_intervals_pos[k1][0]):
                interval_start=significant_GCskew_intervals_pos[k1][0]
                interval_end=min(significant_TAskew_intervals_pos[k2][1],significant_GCskew_intervals_pos[k1][1])
                # if interval_end-interval_start>=interval_threshold:
                overlapping_pos_intervals.append((interval_start,interval_end))
    
    # Finding intervals where GC and TA skews are both of negative sign, and overlap 
    overlapping_neg_intervals=[]
    for k1 in range(0,len(significant_GCskew_intervals_neg)):
        for k2 in range(0,len(significant_TAskew_intervals_neg)):
            # case where the start of TA region is withon the GC region
            if (significant_GCskew_intervals_neg[k1][0]<=significant_TAskew_intervals_neg[k2][0]) and (significant_GCskew_intervals_neg[k1][1]>significant_TAskew_intervals_neg[k2][0]):
                interval_start=significant_TAskew_intervals_neg[k2][0]
                interval_end=min(significant_GCskew_intervals_neg[k1][1],significant_TAskew_intervals_neg[k2][1])
                # if interval_end-interval_start>=interval_threshold:
                overlapping_neg_intervals.append((interval_start,interval_end))
            #case where the start of GC region is withon the TA region
            elif (significant_TAskew_intervals_neg[k2][0]<=significant_GCskew_intervals_neg[k1][0]) and (significant_TAskew_intervals_neg[k2][1]>significant_GCskew_intervals_neg[k1][0]):
                interval_start=significant_GCskew_intervals_neg[k1][0]
                interval_end=min(significant_TAskew_intervals_neg[k2][1],significant_GCskew_intervals_neg[k1][1])
                # if interval_end-interval_start>=interval_threshold:
                overlapping_neg_intervals.append((interval_start,interval_end))

    return overlapping_pos_intervals,overlapping_neg_intervals

def define_GRINS_predicted_intervals_AUC(pos_intervals,neg_intervals,AUC_threshold):
    predicted_GRINS_AUC=[]
    for i in range(0,len(pos_intervals)):
        overlapping_interval=pos_intervals[i]
        interval_length=float(overlapping_interval[1]-overlapping_interval[0])
        x=range(overlapping_interval[0],overlapping_interval[1]+1)
        y_GC=GCskew_array[overlapping_interval[0]:overlapping_interval[1]+1]
        y_TA=TAskew_array[overlapping_interval[0]:overlapping_interval[1]+1]
        AUC_GCskew=abs(metrics.auc(x,y_GC))
        AUC_TAskew=abs(metrics.auc(x,y_TA))
        if AUC_GCskew>=AUC_threshold and AUC_TAskew>=AUC_threshold:
            predicted_GRINS_AUC.append(overlapping_interval)
    for i in range(0,len(neg_intervals)):
        overlapping_interval=neg_intervals[i]
        interval_length=float(overlapping_interval[1]-overlapping_interval[0])
        x=range(overlapping_interval[0],overlapping_interval[1]+1)
        y_GC=GCskew_array[overlapping_interval[0]:overlapping_interval[1]+1]
        y_TA=TAskew_array[overlapping_interval[0]:overlapping_interval[1]+1]
        AUC_GCskew=abs(metrics.auc(x,y_GC))
        AUC_TAskew=abs(metrics.auc(x,y_TA))
        if AUC_GCskew>=AUC_threshold and AUC_TAskew>=AUC_threshold:
            predicted_GRINS_AUC.append(overlapping_interval)
    return predicted_GRINS_AUC


def define_nonGRINS_predicted_intervals(pos_intervals,neg_intervals,predicted_GRINS_sorted):
    predicted_nonGRINS=[]
    for item in pos_intervals:
        if item not in predicted_GRINS_sorted:
            predicted_nonGRINS.append(item)
    for item in neg_intervals:
        if item not in predicted_GRINS_sorted:
            predicted_nonGRINS.append(item)
    return predicted_nonGRINS



if __name__ == "__main__":
    # Making a list of annotated GenBank files

    parser = argparse.ArgumentParser(description='Get the range of cluster IDs to parse')

    parser.add_argument('-startID','--start_ID')
    parser.add_argument('-endID','--end_ID')
    parser.add_argument('-inputf','--input_folder')
    parser.add_argument('-outputf','--output_folder')

    window=150
    step=30
    AUC_threshold=6

    args=parser.parse_args()

    my_in_path=args.input_folder
    my_out_path=args.output_folder

    files = [f for f in listdir(my_in_path) if isfile(join(my_in_path, f))]
    GRINS_positives=0
    # GC_contents=[]


    cluster_ID_start=int(args.start_ID)
    cluster_ID_end=args.end_ID

    if cluster_ID_end=="3552":
        cluster_ID_end=len(files)
    else:
        cluster_ID_end=int(cluster_ID_end)

    current_files=files[cluster_ID_start:cluster_ID_end]

    with open(my_out_path+"Detection_results_orphans"+str(cluster_ID_start)+"_"+str(cluster_ID_end)+".txt",'w') as results_file:
        for item in current_files:
            
            record = SeqIO.read(my_in_path+item, "gb")
            input_DNA=record.seq
            record_name=record.name

            cluster_num=item.split(".")[1]
            if cluster_num[7:9]=="00":
                cluster_n=cluster_num[9]
            elif cluster_num[7]=="0":
                cluster_n=cluster_num[8:]
            else:
                cluster_n=cluster_num[7:]

            # Calculating skews
            GCskew_array=[]
            TAskew_array=[]
            for i in range(0,len(input_DNA),step):
                fragment=input_DNA[i:i+window]
                GCskew_array.append(GC_skew(fragment))
                TAskew_array.append(TA_skew(fragment))

            GC_content=np.around(float(input_DNA.count("G")+input_DNA.count("C"))/float(len(input_DNA)),4)
            # GC_contents.append(GC_content)

            # Finding intervals where GC skew crosses the 0 axis
            pos_intervals,neg_intervals=x_crossing(GCskew_array,TAskew_array)
            if len(pos_intervals)==0 and len(neg_intervals)==0:
                print "Problem occurred cor cluster",record_name,cluster_n
                results_file.write(record_name+" "+cluster_n+"\t"+str(GC_content)+"\t"+"?\n")
            else:

                # Calculating metrics for the AUC threshold

                predicted_GRINS_AUC=define_GRINS_predicted_intervals_AUC(pos_intervals,neg_intervals,AUC_threshold)
                predicted_GRINS_sorted_AUC=sorted(predicted_GRINS_AUC, key=lambda tup: tup[0])

                # predicted_nonGRINS_AUC=define_nonGRINS_predicted_intervals(pos_intervals,neg_intervals,predicted_GRINS_sorted_AUC)
                # predicted_nonGRINS_sorted_AUC=sorted(predicted_nonGRINS_AUC, key=lambda tup: tup[0])

                if len(predicted_GRINS_sorted_AUC)>0:
                    GRINS_positives+=1
                    results_file.write(record_name+" "+cluster_n+"\t"+str(GC_content)+"\t"+"GRINS\n")
                else:
                    results_file.write(record_name+" "+cluster_n+"\t"+str(GC_content)+"\t"+"-\n")

    print "The number of GRINS-positive PKSs is:",GRINS_positives



