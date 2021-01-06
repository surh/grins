#!/usr/bin/env python
# Copyright (C) 2019 Aleksandra Nivina

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# To see the GNU General Public License, please visit
# <http://www.gnu.org/licenses/>.

import argparse
import os
from os import listdir
from os.path import isfile, isdir, join
import numpy as np
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
import matplotlib.colors as clr
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
from Bio.SeqFeature import SeqFeature,FeatureLocation
from scipy.stats import pearsonr


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("This script detects GRINS"
                          "(genetic repeats of intense nulceotide skews)"
                          "based on sequence duplication and skew patterns,"
                          "in a folder containing GenBank (default) or Fasta files.")


    required.add_argument("--gb_input", help=("GenBank Input folder."),
                          required=True, type=str)

    required.add_argument("--startID", help=("GenBank Input folder."),
                          required=True, type=int)

    required.add_argument("--endID", help=("GenBank Input folder."),
                          required=True, type=int)

    required.add_argument("--similarity_input", help=("Similarity Input folder."),
                          required=True, type=str)

    # # Define other arguments
    parser.add_argument("--plot_output", help=("Name of the plot output folder."
                        "If '', then folder 'GRINS_plots'"
                        "is created in the current directory."),
                        type=str,
                        default='./GRINS_plots')
    parser.add_argument("--gb_output", help=("Name of the GenBank output folder."
                    "If '', then folder 'GRINS_results'"
                    "is created in the current directory."),
                    type=str,
                    default='./GRINS_results')
    parser.add_argument("--format", help=("Format of the input files: GenBank or Fasta."),
                        type=str,
                        default='gb',
                        choices=['gb', 'fasta'])
    parser.add_argument("--skew_window", help=("Sliding window size."
                        "Default: 150."),
                        type=int,
                        default=150)
    parser.add_argument("--skew_step", help=("Sliding window step."
                        "Default: 30."),
                        type=int,
                        default=30)
    parser.add_argument("--similarity_threshold", help=("Threshold for the maximum pairwise similarity."
                        "Default: 0.8"),
                        type=float,
                        default=0.8)

    parser.add_argument("--intensity_threshold", help=("Threshold for the mean skew intensity."
                        "Default: 0.15"),
                        type=float,
                        default=0.15)

    args = parser.parse_args()

    return args


# Calculating GC skews
def GC_skew(input_DNA,window,step):
    GCskew_array=[]
    for i in range(0,len(input_DNA)-window,step):
        fragment=input_DNA[i:i+window]
        G_number=fragment.count("G")
        C_number=fragment.count("C")
        if G_number+C_number!=0:
            skew=float(G_number-C_number)/float(G_number+C_number)
        else:
            skew=0
        GCskew_array.append(skew)

    seq_length=len(input_DNA)
    array_length=(seq_length-window)//30
    GCskew_array=GCskew_array[:array_length]

    return GCskew_array


# Calculating TA skews
def TA_skew(input_DNA,window,step):
    TAskew_array=[]
    for i in range(0,len(input_DNA)-window,step):
        fragment=input_DNA[i:i+window]
        T_number=fragment.count("T")
        A_number=fragment.count("A")
        if T_number+A_number!=0:
            skew=float(T_number-A_number)/float(T_number+A_number)
        else:
            skew=0
        TAskew_array.append(skew)

    seq_length=len(input_DNA)
    array_length=(seq_length-window)//30
    TAskew_array=TAskew_array[:array_length]

    return TAskew_array


# Calculating maximum similarity scores
def finding_max_similarity(record_name,seq_length,window,filename,similarity_value_folder):
    with open(similarity_value_folder+filename+".txt", 'r') as similarity_file:
        similarity_score_list=similarity_file.readlines()
    max_score_array=[]
    for i in range(0,len(similarity_score_list)):
        line=similarity_score_list[i].strip("\n").split(",")
        max_score=0
        for j in range(0,len(line)):
            if j<=i-5 or j>=i+5:
                score=float(line[j])/150.0
                if score>max_score:
                    max_score=score
        max_score_array.append(max_score)

    array_length=(seq_length-window)//30
    max_score_array=max_score_array[:array_length]

    return max_score_array


# Finding regions of high sequence similarity
def finding_duplicated_regions_pairwise(similarity_max_array,similarity_threshold,seq_length,window):
    dupl_region_starts=[]
    dupl_region_ends=[]
    
    # this list will contain absolute coordinates of duplicated regions
    dupl_regions=[]
    
    # this array will contain "yes"/"no" value for each window, to represent whether it is duplicated or not
    dupl_region_array=[]

    array_length=(seq_length-window)//30
    # print len(similarity_max_array),array_length

    # finding duplicated regions (window cooredinates) based on sequence similarity
    for i in range(0,array_length):
        if len(dupl_region_starts)==len(dupl_region_ends):
            if np.mean(similarity_max_array[np.max([i-5,0]):i])<similarity_threshold and np.mean(similarity_max_array[i:np.min([i+5,array_length])])>=similarity_threshold:
                dupl_region_starts.append(i)
        else:
            if np.mean(similarity_max_array[i:np.min([i+5,array_length])])<similarity_threshold:
                dupl_region_ends.append(i)

    # adding last region's end, in case it is the end of the sequence
    if len(dupl_region_ends)<len(dupl_region_starts):
        dupl_region_ends.append(array_length)

    # combining overlapping regions and only keeping regions >500bp
    if len(dupl_region_starts)>0:
        prev_start=dupl_region_starts[0]
        prev_end=dupl_region_ends[0]

        for k in range(1,len(dupl_region_starts)):
            start=dupl_region_starts[k]
            end=dupl_region_ends[k]

            # combine previous and current regions if they overlap
            if start*30 <= prev_end*30+150:
                prev_start=prev_start
                prev_end=end

            else:
                # record previous region, if it is longer than 500bp
                # if prev_end*30+150 - prev_start*30 >= 500:
                dupl_regions.append((prev_start,prev_end))

                # updating what is considered to be previous region
                prev_start=start
                prev_end=end


        # saving the last region if it is longer than 500bp
        # if prev_end*30+150 - prev_start*30 >= 500:
        dupl_regions.append((prev_start,prev_end))


    # generating a "yes"/"no" array that tells if each window is part of a duplicated region or not    
    if len(dupl_regions)>0:
        # first negative region
        dupl_region_array.extend(["no"]* dupl_regions[0][0] )
        # first positive region
        dupl_region_array.extend(["yes"]* (dupl_regions[0][1]-dupl_regions[0][0]) )
    
        # the rest
        for k in range(1,len(dupl_regions)):
            dupl_region_array.extend(["no"]* (dupl_regions[k][0]-dupl_regions[k-1][1]) )
            dupl_region_array.extend(["yes"]* (dupl_regions[k][1]-dupl_regions[k][0]) )

        # last negative
        dupl_region_array.extend(["no"]* (array_length-dupl_regions[-1][1]) )


    # case where there are no duplicated regions
    else:
        dupl_region_array.extend(["no"]* array_length)


    if len(dupl_region_array)!=array_length:
        print "boolean array length",len(dupl_region_array),"similarity array",array_length

    return dupl_regions,dupl_region_array



def finding_correlations(GCskew_array,TAskew_array,skew_window,skew_step,dupl_regions,seq_length,window):
    dupl_region_correlations=[]
    correlations_array=[]
    pvalue_array=[]
    array_length=(seq_length-window)//30

    # case where there are any duplicated regions
    if len(dupl_regions)>0:
        prev_region_end=0
        prev_skew_region_end=0

        for i in range(0,len(dupl_regions)):
            region_start=dupl_regions[i][0]
            skew_region_start=region_start*30//skew_step
            region_end=dupl_regions[i][1]
            skew_region_end=region_end*30//skew_step

            # correlation in the negative region
            GC_skew_values=GCskew_array[prev_skew_region_end:skew_region_start]
            TA_skew_values=TAskew_array[prev_skew_region_end:skew_region_start]
            correlation,p_value=pearsonr(GC_skew_values,TA_skew_values)
            if str(correlation)=="nan":
                correlation=0
                p_value=1
            
            # saving values in the array for similarity value windows, not skew value windows
            correlations_array.extend([correlation]* (region_start-prev_region_end) )
            pvalue_array.extend([p_value]* (region_start-prev_region_end) )

            # correlations in positive region
            GC_skew_values=GCskew_array[skew_region_start:skew_region_end]
            TA_skew_values=TAskew_array[skew_region_start:skew_region_end]
            correlation,p_value=pearsonr(GC_skew_values,TA_skew_values)
            if str(correlation)=="nan":
                correlation=0
                p_value=1
            
            # saving the correlation value for this duplicated region
            dupl_region_correlations.append(correlation)

            # saving values in the array for similarity value windows, not skew value windows
            correlations_array.extend([correlation]* (region_end-region_start) )
            pvalue_array.extend([p_value]* (region_end-region_start) )
            prev_region_end=region_end

        # adding the last negative region to the array
        GC_skew_values=GCskew_array[prev_region_end:len(GCskew_array)]
        TA_skew_values=TAskew_array[prev_region_end:len(GCskew_array)]
        correlation,p_value=pearsonr(GC_skew_values,TA_skew_values)
        if str(correlation)=="nan":
            correlation=0
            p_value=1
        
        # saving values in the array for similarity value windows, not skew value windows  
        correlations_array.extend([correlation]* (array_length-prev_region_end) )
        pvalue_array.extend([p_value]* (array_length-prev_region_end) )

    # case where there are no duplicated regions
    else:
        GC_skew_values=GCskew_array
        TA_skew_values=TAskew_array
        correlation,p_value=pearsonr(GC_skew_values,TA_skew_values)
        if str(correlation)=="nan":
            correlation=0
            p_value=1
        
        correlations_array=[correlation]* array_length
        pvalue_array=[p_value]* array_length

    if len(correlations_array)!=array_length:
        print "correlation array lenth",len(correlations_array),"similariry array",array_length

    return dupl_region_correlations,correlations_array,pvalue_array



def finding_mean_skew(GCskew_array,TAskew_array,skew_window,skew_step,dupl_regions,seq_length,window):
    GCskew_mean_array=[]
    TAskew_mean_array=[]
    dupl_GCskew_means=[]
    dupl_TAskew_means=[]
    array_length=(seq_length-window)//30

    # case where there are any duplicated regions
    if len(dupl_regions)>0:
        prev_region_end=0
        prev_skew_region_end=0

        for i in range(0,len(dupl_regions)):
            region_start=dupl_regions[i][0]
            skew_region_start=region_start*30//skew_step
            region_end=dupl_regions[i][1]
            skew_region_end=region_end*30//skew_step

            # mean skews in the negative region
            GC_skew_values=np.abs(GCskew_array[prev_skew_region_end:skew_region_start])
            TA_skew_values=np.abs(TAskew_array[prev_skew_region_end:skew_region_start])
            if len(GC_skew_values)>0:
                GCskew_mean=np.mean(GC_skew_values)
                TAskew_mean=np.mean(TA_skew_values)
            else:
                GCskew_mean=0
                TAskew_mean=0
            
            # saving values in the array for similarity value windows, not skew value windows
            GCskew_mean_array.extend([GCskew_mean]* (region_start-prev_region_end) )
            TAskew_mean_array.extend([TAskew_mean]* (region_start-prev_region_end) )

            # mean skews in positive region
            GC_skew_values=np.abs(GCskew_array[skew_region_start:skew_region_end])
            TA_skew_values=np.abs(TAskew_array[skew_region_start:skew_region_end])
            if len(GC_skew_values)>0:
                GCskew_mean=np.mean(GC_skew_values)
                TAskew_mean=np.mean(TA_skew_values)
            else:
                GCskew_mean=0
                TAskew_mean=0

            # saving values in the array for similarity value windows, not skew value windows
            GCskew_mean_array.extend([GCskew_mean]* (region_end-region_start) )
            TAskew_mean_array.extend([TAskew_mean]* (region_end-region_start) )
            dupl_GCskew_means.append(GCskew_mean)
            dupl_TAskew_means.append(TAskew_mean)
            prev_region_end=region_end

        # adding mean skews of the last negative region to the  array
        GC_skew_values=np.abs(GCskew_array[prev_skew_region_end:len(GCskew_array)])
        TA_skew_values=np.abs(TAskew_array[prev_skew_region_end:len(GCskew_array)])
        if len(GC_skew_values)>0:
            GCskew_mean=np.mean(GC_skew_values)
            TAskew_mean=np.mean(TA_skew_values)
        else:
            GCskew_mean=0
            TAskew_mean=0
        
        # saving values in the array for similarity value windows, not skew value windows
        GCskew_mean_array.extend([GCskew_mean]* (array_length-prev_region_end) )
        TAskew_mean_array.extend([TAskew_mean]* (array_length-prev_region_end) )

    # case where there are no duplicated regions
    else:
        GC_skew_values=np.abs(GCskew_array)
        TA_skew_values=np.abs(TAskew_array)
        if len(GC_skew_values)>0:
            GCskew_mean=np.mean(GC_skew_values)
            TAskew_mean=np.mean(TA_skew_values)
        else:
            GCskew_mean=0
            TAskew_mean=0
        
        # saving values in the array for similarity value windows, not skew value windows
        GCskew_mean_array=[GCskew_mean]* array_length
        TAskew_mean_array=[TAskew_mean]*array_length

    if len(GCskew_mean_array)!=array_length:
        print "Mean GC skew array lenth",len(GCskew_mean_array),"similarity array",array_length

    return dupl_GCskew_means,dupl_TAskew_means,GCskew_mean_array,TAskew_mean_array



# Discarding duplicated regions with low absolute correlations
def finding_GRINS(dupl_regions,pos_GCskew_means,pos_TAskew_means,seq_length,intensity_threshold,window):
    GRINS_regions=[]
    GRINS_regions_array=[]

    array_length=(seq_length-window)//30

    # print len(dupl_regions),len(pos_GCskew_means),len(pos_TAskew_means)
    for i in range(0,len(pos_GCskew_means)):
        if pos_GCskew_means[i]>=intensity_threshold and pos_TAskew_means[i]>=intensity_threshold and dupl_regions[i][1]*30+150-dupl_regions[i][0]*30>500:
            GRINS_regions.append(dupl_regions[i])

    if len(GRINS_regions)==0:
        GRINS_regions_array=["no"]*array_length
    else:

        if GRINS_regions[0][0]==0:
            GRINS_regions_array.extend(["yes"]*(GRINS_regions[0][1]-GRINS_regions[0][0]))
        else:
            GRINS_regions_array.extend(["no"]*(GRINS_regions[0][0]))
            GRINS_regions_array.extend(["yes"]*(GRINS_regions[0][1]-GRINS_regions[0][0]))

        for n in range(1,len(GRINS_regions)):
            GRINS_regions_array.extend(["no"]*(GRINS_regions[n][0]-GRINS_regions[n-1][1]))
            GRINS_regions_array.extend(["yes"]*(GRINS_regions[n][1]-GRINS_regions[n][0]))

        if GRINS_regions[-1][1]!=array_length:
            GRINS_regions_array.extend(["no"]*(array_length-GRINS_regions[-1][1]))

    if len(GRINS_regions_array)!=array_length:
        print "Mean GC skew array lenth",len(GRINS_regions_array),"skew array length",array_length

    return GRINS_regions, GRINS_regions_array



# Plotting the graph for each record
def plot_graph(record_name,plot_folder,max_similarity_array,GCskew_array,TAskew_array,dupl_regions,GRINS_regions,window,step,seq_length):
        a=np.arange(0,len(max_similarity_array)*30,30)
        plt.figure(0)
        fig,ax1=plt.subplots(figsize=(10, 4))

        ax1.set_ylim(0,1)
        ax1.set_xlim(0,seq_length)
        ax1.grid(False)
        plt.plot(a,max_similarity_array,color="black",linewidth=0.5)
        plt.legend(["Pairwise similarity"],loc=3)
        for j in range(0,len(dupl_regions)):
            region_start=dupl_regions[j][0]*30
            region_end=dupl_regions[j][1]*30
            rect = patches.Rectangle((region_start,1),region_end-region_start,-2,linewidth=1,edgecolor='grey',facecolor='grey', alpha=0.3)
            ax1.add_patch(rect)
        for k in range(0,len(GRINS_regions)):
            region_start=GRINS_regions[k][0]*30
            region_end=GRINS_regions[k][1]*30
            rect = patches.Rectangle((region_start,1),region_end-region_start,-2,linewidth=1,edgecolor='teal',facecolor='teal', alpha=0.3)
            ax1.add_patch(rect)

        ax2=ax1.twinx()
        ax2.set_ylim(-1,1)
        ax2.set_xlim(0,seq_length)
        ax2.grid(False)
        plt.plot(a,GCskew_array,color="blue",linewidth=0.5)
        plt.plot(a,TAskew_array,color="red",linewidth=0.5)
        plt.legend(["GC skew","TA skew"],loc=4)

        plt.title("GRINS in %s" % (record_name) )
        plt.savefig("./%s/%s_pairwise_similarity_%s.png" % (plot_folder,record_name,similarity_threshold) )
        plt.close('all')



# Save GenBank file with duplicated and GRINS regions annotated
def save_annotated_Genbank(record,genbank_output_folder,genbank_output_name,dupl_regions,GRINS_regions):
    with open("%s/%s.gbk" %(genbank_output_folder,genbank_output_name),"w") as annotated_file:
        record_annotated=record
        for i in range(0,len(dupl_regions)):
            feature=SeqFeature(FeatureLocation(start=dupl_regions[i][0]*30, end=dupl_regions[i][1]*30+150), type='Duplicated')
            record_annotated.features.append(feature)
        for i in range(0,len(GRINS_regions)):
            feature=SeqFeature(FeatureLocation(start=GRINS_regions[i][0]*30, end=GRINS_regions[i][1]*30+150), type='GRINS')
            record_annotated.features.append(feature)           
        SeqIO.write(record_annotated, annotated_file, 'genbank')


if __name__ == "__main__":


    # Reading arguments and creating result folders
    args = process_arguments()

    mypath = args.gb_input
    similarity_value_folder=args.similarity_input

    gb_folder= args.gb_output

    plot_folder=args.plot_output

    cluster_ID_start=args.startID
    cluster_ID_end=args.endID

    skew_window=args.skew_window
    skew_step=args.skew_step
    similarity_threshold=args.similarity_threshold
    intensity_threshold=args.intensity_threshold

    files = [f for f in listdir(mypath) if isfile(join(mypath, f))]
    # files.remove("CP025990.cluster010.gbk")
    # files.remove("CP025990.cluster011.gbk")

    # print cluster_ids
    if cluster_ID_end==3551:
        cluster_ID_end=len(files)
    else:
        cluster_ID_end=cluster_ID_end
    # print len(cluster_ids)
    cluster_ids_curr=files[cluster_ID_start:cluster_ID_end]
    # print "Started parsing clusters",cluster_ids_curr


    # Opening results file
    with open("./GRINS_detection_results_pairwise%s_w%s_s%s_intensity%s_%s-%s.txt" %(similarity_threshold,skew_window,skew_step,intensity_threshold,cluster_ID_start,cluster_ID_end),'w') as results_file:
        with open("./Values_pairwise%s_w%s_s%s_intensity%s.txt" %(similarity_threshold,skew_window,skew_step,intensity_threshold),'w') as values_output_file:
            values_output_file.write("\t".join(["Record","GC skew","TA skew","Skew correlation","Correlation p-value","Mean GC skew","Mean TA skew","Max similarity","Duplicated","GRINS"])+"\n")

        
            with open("List_of_missing_clusterfiles_%s-%s.txt" %(cluster_ID_start,cluster_ID_end),'w') as missing_cluster_file:
                for filename in cluster_ids_curr:

                    filename_no_extension=filename[:filename.find(".gb")]
                
                    try:
                        open(similarity_value_folder+filename_no_extension+".txt",'r')
                    except:
                        missing_cluster_file.write(filename_no_extension+".txt\n")
                    else:

                        record = SeqIO.read(mypath+filename, "gb")        

                        input_DNA=record.seq
                        input_DNA_str=str(input_DNA)
                        record_name=record.name
                        cluster_name=filename
                        print "Working on record",record.id


                        # Calculating skews
                        GCskew_array=GC_skew(input_DNA,skew_window,skew_step)
                        TAskew_array=TA_skew(input_DNA,skew_window,skew_step)

                        # Reading the pairwise alignment file and finding maximum similarity score
                        max_similarity_array=finding_max_similarity(record_name,len(input_DNA_str),skew_window,filename_no_extension,similarity_value_folder)

                        # STEP1: Finding duplicated regions (maximum pairwise similarity > similarity_threshold)
                        dupl_regions,dupl_regions_array=finding_duplicated_regions_pairwise(max_similarity_array,similarity_threshold,len(input_DNA_str),skew_window)

                        # Calculating correlations between GC and TA skews for each region
                        pos_region_correlations,correlations_array,pvalue_array=finding_correlations(GCskew_array,TAskew_array,skew_window,skew_step,dupl_regions,len(input_DNA_str),skew_window)

                        # Calculating mean absolute GC and TA skews for each region
                        pos_GCskew_means,pos_TAskew_means,GCskew_mean_array,TAskew_mean_array=finding_mean_skew(GCskew_array,TAskew_array,skew_window,skew_step,dupl_regions,len(input_DNA_str),skew_window)

                        # STEP 2: Finding GRINS (mean skew intensity > intensity_threshold)
                        GRINS_regions, GRINS_regions_array=finding_GRINS(dupl_regions,pos_GCskew_means,pos_TAskew_means,len(input_DNA_str),intensity_threshold,skew_window)


                        # Writing result into file, if GRINS found
                        results_file.write(cluster_name+"\t")
                        if len(dupl_regions)>0:
                            results_file.write("Duplications found:\t")
                            results_file.write("; ".join([str(k[0]*30)+"-"+str(k[1]*30+150) for k in dupl_regions])+"\t")
                        else:
                            results_file.write("Duplications not found.\t\t")
                        if len(GRINS_regions)>0:
                            results_file.write("GRINS found:\t")
                            results_file.write("; ".join([str(k[0]*30)+"-"+str(k[1]*30+150) for k in GRINS_regions])+"\n")
                        else:
                            results_file.write("GRINS not found.\t\n")

                        # Saving the annotated GenBank file
                        genbank_output_name=filename+".gb"
                        save_annotated_Genbank(record,gb_folder,filename_no_extension,dupl_regions,GRINS_regions)

                        # Plotting the skews and the number of duplicated k-mers
                        plot_graph(cluster_name,plot_folder,max_similarity_array,GCskew_array,TAskew_array,dupl_regions,GRINS_regions,skew_window,skew_step,len(input_DNA_str))

                        # # Outputting different values for each window, to make a pairs plot later
                        for i in range(0,len(max_similarity_array)):
                            # print len(GCskew_array),len(max_similarity_array),len(correlations_array),len(dupl_regions_array),len(GRINS_regions_array)
                            values_output_file.write("\t".join([cluster_name,str(GCskew_array[i]),str(TAskew_array[i]),str(correlations_array[i]),str(pvalue_array[i]),str(GCskew_mean_array[i]),str(TAskew_mean_array[i]),str(max_similarity_array[i]),str(dupl_regions_array[i]),str(GRINS_regions_array[i])])+"\n")
                            if max_similarity_array[i]<0.6:
                                print "Low similarity values",cluster_name, i



