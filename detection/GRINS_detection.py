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
from sklearn import metrics


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

    # Define required arguments
    required.add_argument("--input", help=("Input folder."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Name of the output folder."
                        "If '', then folder 'GRINS_results'"
                        "is created in the current directory."),
                        type=str,
                        default='./GRINS_results')
    parser.add_argument("--format", help=("Format of the input files: GenBank or Fasta."),
                        type=str,
                        default='gb',
                        choices=['gb', 'fasta'])
    parser.add_argument("--window", help=("Sliding window size."
                        "Default: 150."),
                        type=int,
                        default=150)
    parser.add_argument("--step", help=("Sliding window step."
                        "Default: 30."),
                        type=int,
                        default=30)
    parser.add_argument("--kmer", help=("kmer size for detecting duplications."
                        "Default: 11."),
                        type=int,
                        default=11)
    parser.add_argument("--kmer_threshold", help=("Threshold for the number of duplicated kmers in the window."
                        "Default: 70."),
                        type=int,
                        default=70)
    parser.add_argument("--skew_threshold", help=("Threshold for the AUC under the curves of GC and TA skews."
                        "Default: 0.5."),
                        type=float,
                        default=0.5)

    args = parser.parse_args()

    return args


# Calculating GC skews
def GC_skew(input_DNA,window,step):
    GCskew_array=[]
    for i in range(0,len(input_DNA)-window+1,step):
        fragment=input_DNA[i:i+window]
        G_number=fragment.count("G")
        C_number=fragment.count("C")
        if G_number+C_number!=0:
            skew=float(G_number-C_number)/float(G_number+C_number)
        else:
            skew=0
        GCskew_array.append(skew)

    return GCskew_array


# Calculating TA skews
def TA_skew(input_DNA,window,step):
    TAskew_array=[]
    for i in range(0,len(input_DNA)-window+1,step):
        fragment=input_DNA[i:i+window]
        T_number=fragment.count("T")
        A_number=fragment.count("A")
        if T_number+A_number!=0:
            skew=float(T_number-A_number)/float(T_number+A_number)
        else:
            skew=0
        TAskew_array.append(skew)

    return TAskew_array


# Generating existing kmers for each window in the input sequence
def generating_kmer_dictionary(seq,start,step,window,k):
    Kmers=dict()
    if type(seq) is not str:
        raise ValueError("seq must be of type 'str'.")
    if step < 1 or window < 1 or k < 1:
        raise ValueError("step, window and k must be positive integers.")


    # going through the sequence in forward direction
    for w_start in range(start, len(seq) - k + 1):
        kmer = seq[w_start:w_start+k]

        if kmer in Kmers:
            Kmers[kmer] = Kmers[kmer] + 1
        else:
            Kmers[kmer] = 1
    
    # making reverse complement
    sequence=Seq(seq)
    sequence_revcomp=sequence.reverse_complement()
    seq_revcomp=str(sequence_revcomp)

    # going through the sequence in reverse direction
    for w_start in range(0, len(seq) - k - start + 1):
        kmer = seq_revcomp[w_start:w_start+k]

        if kmer in Kmers:
            Kmers[kmer] = Kmers[kmer] + 1
        else:
            Kmers[kmer] = 1

    return Kmers


# Calculating the number of duplicated kmers for each window in the input sequence
def duplicated_kmers(seq,start,step,window,k,Kmers):

    if type(seq) is not str:
        raise ValueError("seq must be of type 'str'.")
    if step < 1 or window < 1 or k < 1:
        raise ValueError("step, window and k must be positive integers.")

    
    duplications_seq = []
    # calculating the number of duplicated kmers in each window, using kmers found in the entire DNA sequence, on both strands
    for w_start in range(start, len(seq) - window + 1, step):
        w_end = min(w_start + window, len(seq))

        curr_kmers=[]

        for start in range(w_start, w_end - k + 1):
            kmer = seq[start:start+k]
            if kmer not in curr_kmers:
                curr_kmers.append(kmer)

        duplicated_kmers=np.sum(1 for x in curr_kmers if Kmers[x]>1)
        duplications_seq.append(duplicated_kmers)

    return duplications_seq

# Detecting duplicated regions, based on kmers
def duplicated_regions(duplicated_kmer_array,kmer_threshold):
    dupl_region_starts=[]
    dupl_region_ends=[]
    for i in range(0,len(duplicated_kmer_array)):

        # if currently not in duplicated region
        # and sequence similarity in current position and in next 10 positions (averaged) is above threshold:
        # mark this is the start of a duplicated region
        if len(dupl_region_starts)==len(dupl_region_ends):
            if duplicated_kmer_array[i]>=kmer_threshold and np.mean(duplicated_kmer_array[i:i+5])>=kmer_threshold:
                dupl_region_starts.append(i)

        # if currently in duplicated region
        # and sequence similarity in current position is above threshold, but sequence similarity in next 10 positions (averaged) is above threshold:
        # mark this is the end of a duplicated region
        else:
            if duplicated_kmer_array[i]>=kmer_threshold and np.mean(duplicated_kmer_array[i+1:i+6])<kmer_threshold:
                dupl_region_ends.append(i)

    # if no duplicated region end found so far, consider the end of the sequence as the end of the duplicated region
    if len(dupl_region_starts)!=len(dupl_region_ends):
        dupl_region_ends.append(len(duplicated_kmer_array))

    return dupl_region_starts,dupl_region_ends


# Calculating absolute AUCs under GC skew and TA skew plots 
def duplicated_region_skews(dupl_region_starts,dupl_region_ends,GCskew_array,TAskew_array):
    AUC_GCskew_regions=[]
    AUC_TAskew_regions=[]
    for r in range(0,len(dupl_region_starts)):

        # if duplicated regions includes more than 1 window, calculate AUC of skew function
        if dupl_region_ends[r]-dupl_region_starts[r]>1:
            region_GCskew=np.abs(GCskew_array[dupl_region_starts[r]:dupl_region_ends[r]])
            AUC_GCskew=metrics.auc(np.arange(dupl_region_starts[r],dupl_region_ends[r]),region_GCskew)
            AUC_GCskew_regions.append(AUC_GCskew)

            region_TAskew=np.abs(TAskew_array[dupl_region_starts[r]:dupl_region_ends[r]])
            AUC_TAskew=metrics.auc(np.arange(dupl_region_starts[r],dupl_region_ends[r]),region_GCskew)
            AUC_TAskew_regions.append(AUC_TAskew)

        # otherwise, report skews in one window
        else:
            region_GCskew=np.abs(GCskew_array[dupl_region_starts[r]:dupl_region_ends[r]])
            region_TAskew=np.abs(TAskew_array[dupl_region_starts[r]:dupl_region_ends[r]])
            AUC_GCskew_regions.append(region_GCskew[0])
            AUC_TAskew_regions.append(region_TAskew[0])

    return AUC_GCskew_regions,AUC_TAskew_regions


# Among duplicated regions, searching those that have significant skews
def skeweded_duplicated_regions(dupl_region_starts,dupl_region_ends,AUC_GCskew_regions,AUC_TAskew_regions,skew_threshold):
    skewed_dupl_regions=[]
    for r in range(0,len(dupl_region_starts)):
        if (AUC_GCskew_regions[r]>=skew_threshold) and (AUC_TAskew_regions[r]>=skew_threshold): # BOTH skews are intense
        # if (AUC_GCskew_regions[r]>=skew_threshold) or (AUC_TAskew_regions[r]>=skew_threshold): # EITHER of the skews is intense
            skewed_dupl_regions.append([dupl_region_starts[r],dupl_region_ends[r]])

    return skewed_dupl_regions


# Plotting the graph for each record
def plot_graph(record_name,plot_output_folder,duplicated_kmer_array,GCskew_array,TAskew_array,skewed_dupl_regions,window):
        a=np.arange(0,len(duplicated_kmer_array)*step,step)
        plt.figure(0)
        fig,ax1=plt.subplots(figsize=(10, 4))

        ax1.set_ylim(0,window)
        ax1.set_xlim(0,a[-1])
        ax1.grid(False)
        plt.plot(a,duplicated_kmer_array,color="black",linewidth=0.5)
        plt.legend(["Duplicated kmers, k=%d" %kmer],loc=3)
        for j in range(0,len(skewed_dupl_regions)):
            region_start=skewed_dupl_regions[j][0]*step
            region_end=skewed_dupl_regions[j][1]*step
            rect = patches.Rectangle((region_start,window),region_end-region_start,-window,linewidth=1,edgecolor='teal',facecolor='teal', alpha=0.3)
            ax1.add_patch(rect)

        ax2=ax1.twinx()
        ax2.set_ylim(-1,1)
        ax2.set_xlim(0,a[-1])
        ax2.grid(False)
        plt.plot(a,GCskew_array,color="blue",linewidth=0.5)
        plt.plot(a,TAskew_array,color="red",linewidth=0.5)
        plt.legend(["GC skew","TA skew"],loc=4)

        plt.title("GRINS in %s" % (record_name) )
        plt.savefig("%s/%s_GRINS.png" % (plot_output_folder,record_name) )
        plt.close('all')


# Given regions predicted to be GRINS (positive), define other regions of the sequence as negative
def mapping_negative_regions(duplicated_kmer_array,skewed_dupl_regions):
    negative_regions=[]

    # If no duplicated regions, then everything is non-duplicated
    if len(skewed_dupl_regions)==0:
        negative_regions.append([0,len(duplicated_kmer_array)])
    else:

        # Unless duplicated region starts at the very beginning of the sequence, add a first non-duplicated region in the beginning
        if skewed_dupl_regions[0][0]!=0:
            negative_regions.append([0,skewed_dupl_regions[0][0]])

        # If more than one duplicated region, add a non-duplicated region before each of them
        if len(skewed_dupl_regions)>1:
            for m in range(1,len(skewed_dupl_regions)):
                negative_regions.append([skewed_dupl_regions[m-1][1],skewed_dupl_regions[m][0]])

        # Unless duplicated region ends at the very end of the sequence, add a last non-duplicated region in the end
        if skewed_dupl_regions[-1][1]!=len(duplicated_kmer_array):
            negative_regions.append([skewed_dupl_regions[-1][1],len(duplicated_kmer_array)])

    return negative_regions


# Save GenBank file with GRINS regions annotated
def save_annotated_Genbank(record,genbank_output_folder,skewed_dupl_regions):
    with open("%s/%s_GRINS.gb" %(genbank_output_folder,record.name),"w") as annotated_file:
        record_annotated=record
        for i in range(0,len(skewed_dupl_regions)):
            feature=SeqFeature(FeatureLocation(start=skewed_dupl_regions[i][0]*30, end=skewed_dupl_regions[i][1]*30), type='GRINS')
            record_annotated.features.append(feature)
        SeqIO.write(record_annotated, annotated_file, 'genbank')


##############################################################################################
##############################################################################################


if __name__ == "__main__":


    # Reading arguments and creating result folders
    args = process_arguments()

    mypath = args.input

    output_folder= args.output
    if os.path.isdir(args.output):
        raise ValueError("{} file already exists.".format(args.output))
    os.mkdir(output_folder)

    plot_output_folder="%s/Plots" %output_folder
    os.mkdir(plot_output_folder)

    if args.format=="gb":
        genbank_output_folder="%s/Annotated_Genbank_files" %output_folder
        os.mkdir(genbank_output_folder)

    window=args.window
    step=args.step
    kmer=args.kmer
    kmer_threshold=args.kmer_threshold
    skew_threshold=args.skew_threshold
    

    # Opening results file
    with open("%s/GRINS_detection_results.txt" %output_folder,'w') as results_file:
        files = [f for f in listdir(mypath) if isfile(join(mypath, f))]

        # Going over input files
        for item in files:


            # If multiple records per file       
            record_n=len(list(SeqIO.parse(mypath+item, args.format)))
            results_file.write("Out of %d records tested in file %s, GRINS found in the following ones:\n" %(record_n,item))

            
            # The rest is the same, regardless of the data layout
            for record in SeqIO.parse(mypath+item, args.format):
                input_DNA=record.seq
                input_DNA_str=str(input_DNA)
                record_name=record.name
                print "Working on record",record.id

                # Calculating skews
                GCskew_array=GC_skew(input_DNA,window,step)
                TAskew_array=TA_skew(input_DNA,window,step)

                
                # Generating a dictionary of kmers for each record
                Kmers=generating_kmer_dictionary(input_DNA_str,0,step,window,kmer)

                
                # Calculating the number of duplicated kmers per window
                duplicated_kmer_array=duplicated_kmers(input_DNA_str,0,step,window,kmer,Kmers)


                # Detecting duplicated regions, for different kmer duplication number thresholds
                dupl_region_starts,dupl_region_ends=duplicated_regions(duplicated_kmer_array,kmer_threshold)

                
                # Calculating AUC under GC and TA skew graphs (absolute values) for each region
                AUC_GCskew_regions,AUC_TAskew_regions=duplicated_region_skews(dupl_region_starts,dupl_region_ends,GCskew_array,TAskew_array)

                    
                # Among duplicated regions, detecting those that have intense skews, for different skew AUC thresholds
                skewed_dupl_regions=skeweded_duplicated_regions(dupl_region_starts,dupl_region_ends,AUC_GCskew_regions,AUC_TAskew_regions,skew_threshold)


                # Defining negative regions 
                negative_regions=mapping_negative_regions(duplicated_kmer_array,skewed_dupl_regions)

                
                # Writing result into file, if GRINS found
                if len(skewed_dupl_regions)>0:
                    results_file.write("%s\n" %record_name)

                # Plotting the skews and the number of duplicated k-mers
                plot_graph(record_name,plot_output_folder,duplicated_kmer_array,GCskew_array,TAskew_array,skewed_dupl_regions,window)


                # If GenBank files were provided, saving annotated versions 
                if args.format=="gb":
                    save_annotated_Genbank(record,genbank_output_folder,skewed_dupl_regions)



