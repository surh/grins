# import Bio
from Bio import SeqIO
# from Bio.Seq import Seq
# from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature, FeatureLocation
import numpy as np
import argparse
import os
from os import listdir
from os.path import isfile, isdir, join
# import matplotlib
# matplotlib.use('TkAgg')
import matplotlib.pyplot as plt
import matplotlib.patches as patches
# import matplotlib.colors as clr

def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")


    # Define required arguments
    required.add_argument("--seq_input", help=("Folder with GenBank sequences"),
                          required=True, type=str)

    required.add_argument("--dupl_input", help=("Folder with duplication detection results"),
                          required=True, type=str)

    # required.add_argument("--duplBGC_input", help=("Folder with BGC duplication detection results"),
    #                       required=True, type=str)

    parser.add_argument("--seq_output", help=("Annotated GenBank output folder"),
                          required=False, type=str, default="./output/genomes_GRINS")

    parser.add_argument("--GRINS_output", help=("GRINS output folder"),
                          required=False, type=str, default="./output/GRINS.gff3")

    parser.add_argument("--GRINS_BGC_output", help=("GRINS in BGC output folder"),
                          required=False, type=str, default="./output/GRINS_BGC.gff3")

    parser.add_argument("--with_plots", help=("Do you wish to plot GRINS regions? [yes/no]"),
                          required=False, type=str, default="no")

    parser.add_argument("--plot_output", help=("GRINS in BGC output folder"),
                          required=False, type=str, default="./output/plots")
    args = parser.parse_args()

    return args

# Calculating GC skews
def GC_skew(input_DNA,window=150,step=30):
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

    return GCskew_array


# Calculating TA skews
def TA_skew(input_DNA,window=150,step=30):
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

    return TAskew_array

# Calculating mean absolute skews
def abs_skew_means(record,start,end,window=150,step=30):
    input_DNA=record.seq[start:end]
    GCskew_array=[]
    TAskew_array=[]
    for i in range(0,len(input_DNA)-window,step):
        fragment=input_DNA[i:i+window]
        G_number=fragment.count("G")
        C_number=fragment.count("C")
        T_number=fragment.count("T")
        A_number=fragment.count("A")

        if G_number+C_number!=0:
            GCskew=np.abs(float(G_number-C_number)/float(G_number+C_number))
        else:
            GCskew=0
        GCskew_array.append(GCskew)

        if T_number+A_number!=0:
            TAskew=np.abs(float(T_number-A_number)/float(T_number+A_number))
        else:
            TAskew=0
        TAskew_array.append(TAskew)

    mean_GC_skew=np.mean(GCskew_array)
    mean_TA_skew=np.mean(TAskew_array)

    return mean_GC_skew,mean_TA_skew

# Plotting the graph for each record
def plot_graph(plot_folder,assembly,accession,start,end,GCskew_array,TAskew_array,dupl_region_starts,dupl_region_ends,GRINS_starts,GRINS_ends):
        a=np.arange(start,start+len(GCskew_array)*30,30)
        plt.figure(0)
        fig,ax1=plt.subplots(figsize=(10, 4))
        ax1.set_ylim(-1,1)
        ax1.set_xlim(start,end)
        ax1.grid(False)
        plt.plot(a,GCskew_array,color="blue",linewidth=0.5)
        plt.plot(a,TAskew_array,color="red",linewidth=0.5)
        plt.legend(["GC skew","TA skew"],loc=4)
        plt.legend(["Pairwise similarity"],loc=3)
        for j in range(0,len(dupl_region_starts)):
            region_start=dupl_region_starts[j]
            region_end=dupl_region_ends[j]
            rect = patches.Rectangle((region_start,1),region_end-region_start,-2,linewidth=1,edgecolor='grey',facecolor='grey', alpha=0.3)
            ax1.add_patch(rect)
        for k in range(0,len(GRINS_starts)):
            region_start=GRINS_starts[k]
            region_end=GRINS_ends[k]
            rect = patches.Rectangle((region_start,1),region_end-region_start,-2,linewidth=1,edgecolor='teal',facecolor='teal', alpha=0.3)
            ax1.add_patch(rect)


        plt.title("GRINS in %s, %d-%d kb" % (accession,start/1000,end/1000))
        plt.savefig("%s/%s/%s_GRINS_%d-%dkb.png" % (plot_folder,assembly,accession,start/1000,end/1000) )
        plt.close('all')

if __name__ == "__main__":

    args = process_arguments()
    # assembly_accessions=[f[:f.find("_genomic.duplicated.gff3")] for f in listdir(args.dupl_input) if isfile(join(args.dupl_input, f))]
    assembly_accessions=[f[:f.find(".duplicated.gff3")] for f in listdir(args.dupl_input) if isfile(join(args.dupl_input, f))]
    sequence_files=[f for f in listdir(args.seq_input) if isfile(join(args.seq_input, f))]

    if not isdir(args.seq_output):
        os.mkdir(args.seq_output)
    if not isdir(args.GRINS_output):
        os.mkdir(args.GRINS_output)
    if not isdir(args.GRINS_BGC_output):
        os.mkdir(args.GRINS_BGC_output)

    if args.with_plots=="yes":
        if not isdir(args.plot_output):
            os.mkdir(args.plot_output)

    with open("GRINS_detected_in_genomes_and_BGCs.txt",'w') as output_file:
        output_file.write("\t".join(["Genome","Genome length","Number of contigs in the genome","Number of duplications in the genome","Number of BGCSs in the genome","Total length of BGCs in the genome","Total length of T1 PKS BGCs in the genome","Total length of NRPS BGCs in the genome","Total length of terpene BGCs in the genome","Total length of lanthipeptide BGCs in the genome","Total length of ladderane BGCs in the genome","Total length of nucleoside BGCs in the genome","Number of GRINS detected in the genome","Number of GRINS detected in CDSs","Number of GRINS detected in BGCs","Number of GRINS detected in T1 PKS", "Number of GRINS detected in NRPS","Number of GRINS detected in terpene clusters", "Number of GRINS detected in lanthipeptide clusters","Number of GRINS detected in ladderane clusters","Number of GRINS detected in nucleoside clusters"])+"\n")

        for assembly in assembly_accessions:

            length_genome=0
            n_contigs_genome=0

            n_BGC=0
            length_BGC=0

            n_dups=0

            length_PKS=0
            length_NRPS=0
            length_terpene=0
            length_lanthi=0
            length_ladderane=0
            length_nucleoside=0


            n_GRINS_total=0
            n_GRINS_CDS=0
            n_GRINS_BGC=0

            n_GRINS_PKS=0
            n_GRINS_NRPS=0
            n_GRINS_terpenes=0
            n_GRINS_lanthi=0
            n_GRINS_ladderane=0
            n_GRINS_nucleoside=0

            # Reading duplication detection results
            # and creating a dictionary to store these
            # with open("%s/%s_genomic.duplicated.gff3" %(args.dupl_input,assembly),'r') as dupl_results_file:
            with open("%s/%s.duplicated.gff3" %(args.dupl_input,assembly),'r') as dupl_results_file:
                dupl_results=dupl_results_file.readlines()
            dupl_location_dict={}
            for i in range(0,len(dupl_results)):
                if dupl_results[i]!="" and dupl_results[i][0]!="#":
                    n_dups+=1
                    line=dupl_results[i].strip("\n").split("\t")
                    accession=line[0]
                    if accession.find(".")!=-1:
                        accession=accession[:accession.find(".")]
                    start=int(line[3])
                    end=int(line[4])
                    if accession in dupl_location_dict:
                        dupl_location_dict[accession].append([start,end])
                    else:
                        dupl_location_dict[accession]=[[start,end]]

            # Reading sequence file
            filename = [s for s in sequence_files if assembly in s]
            records=SeqIO.parse("./%s/%s" %(args.seq_input,filename[0]),'gb')
            length_var=300

            # Writing a new record, with detected GRINS annotated
            with open("./%s/%s.GRINS.gbk" %(args.seq_output,assembly),"w") as annotated_file:

                # Making several output files for each record, to store some information
                with open("./%s/%s.GRINS_BGC.gff3" %(args.GRINS_BGC_output,assembly),'w') as GRINS_BGC_output_file:
                    GRINS_BGC_output_file.write("\t".join(["Record","GRINS start","GRINS end","BGC"])+"\n")

                    with open("./%s/%s.GRINS_CDS.gff3" %(args.GRINS_output,assembly),'w') as GRINS_CDS_output_file:
                        GRINS_CDS_output_file.write("\t".join(["Record","GRINS start","GRINS end","Locus","CDS name"])+"\n")

                        with open("./%s/%s.GRINS.gff3" %(args.GRINS_output,assembly),'w') as GRINS_output_file:
                            GRINS_output_file.write("\t".join(["Record","GRINS start","GRINS end"])+"\n")


                            # going over records in a GenBank file
                            for record in records:

                                length_genome+=len(record.seq)
                                n_contigs_genome+=1

                                record_features=record.features
                                record_annotated=record
                                accession=record.id
                                if accession.find(".")!=-1:
                                    accession=accession[:accession.find(".")]

                                # counting the numbers and lengths of various BGCs in the genome
                                for feature0 in record_features:
                                    if feature0.type=="region":
                                        BGC_start=np.min([int(feature0.location.start),int(feature0.location.end)])
                                        BGC_end=np.max([int(feature0.location.start),int(feature0.location.end)])
                                        BGC_len=BGC_end-BGC_start
                                        BGC_type=feature0.qualifiers.get("product")[0]

                                        n_BGC+=1
                                        length_BGC+=BGC_len

                                        if BGC_type=="T1PKS" or BGC_type=="transAT-PKS":
                                            length_PKS+=BGC_len
                                        elif BGC_type=="NRPS":
                                            length_NRPS+=BGC_len
                                        elif BGC_type=="terpene":
                                            length_terpene+=BGC_len
                                        elif BGC_type=="lanthipeptide":
                                            length_lanthi+=BGC_len
                                        elif BGC_type=="ladderane":
                                            length_ladderane+=BGC_len
                                        elif BGC_type=="nucleoside":
                                            length_nucleoside+=BGC_len


                                # getting the list of duplicated region from the dictionary created earlier
                                if accession in dupl_location_dict:
                                    dupl_locations=dupl_location_dict[accession]
                                else:
                                    dupl_locations=[]


                                # detecting GRINS using a mean skew threshold
                                GRINS_locations=[]
                                for i in range(0,len(dupl_locations)):

                                    # recording the duplicated region as a feature
                                    dupl_feature=SeqFeature(FeatureLocation(start=dupl_locations[i][0], end=dupl_locations[i][1]), type='Duplication')
                                    record_annotated.features.append(dupl_feature)

                                    if dupl_locations[i][1]-dupl_locations[i][0]>=500:
                                        mean_GCskew,mean_TAskew=abs_skew_means(record,dupl_locations[i][0],dupl_locations[i][1],150,30)

                                        # GRINS found
                                        if mean_GCskew>=0.15 and mean_TAskew>=0.15:

                                            # counting this GRINS and recording it as feature
                                            n_GRINS_total+=1
                                            GRINS_feature=SeqFeature(FeatureLocation(start=dupl_locations[i][0], end=dupl_locations[i][1]), type='GRINS')
                                            record_annotated.features.append(GRINS_feature)
                                            GRINS_start=int(dupl_locations[i][0])
                                            GRINS_end=int(dupl_locations[i][1])
                                            GRINS_locations.append((GRINS_start,GRINS_end))


                                            # checking if GRINS is in a BGC
                                            in_BGC=0
                                            for feature1 in record_features:
                                                if feature1.type=="region":
                                                    curr_BGC_start=np.min([int(feature1.location.start),int(feature1.location.end)])
                                                    curr_BGC_end=np.max([int(feature1.location.start),int(feature1.location.end)])
                                                    curr_BGC_type=feature1.qualifiers.get("product")[0]
                                                    if GRINS_start>=(curr_BGC_start-length_var) and GRINS_start<(curr_BGC_end+length_var) and GRINS_end>(curr_BGC_start-length_var) and GRINS_end<=(curr_BGC_end+length_var):
                                                        # GRINS in BGC
                                                        in_BGC=1
                                                        n_GRINS_BGC+=1

                                                        # checking in what type of BGC this GRINS is in
                                                        if curr_BGC_type=="T1PKS" or curr_BGC_type=="transAT-PKS":
                                                            n_GRINS_PKS+=1
                                                        elif curr_BGC_type=="NRPS":
                                                            n_GRINS_NRPS+=1
                                                        elif curr_BGC_type=="terpene":
                                                            n_GRINS_terpenes+=1
                                                        elif curr_BGC_type=="lanthipeptide":
                                                            n_GRINS_lanthi+=1
                                                        elif curr_BGC_type=="ladderane":
                                                            n_GRINS_ladderane+=1
                                                        elif curr_BGC_type=="nucleoside":
                                                            n_GRINS_nucleoside+=1
                                                        break

                                            # checking if GRINS is in a CDS
                                            in_CDS=0
                                            for feature2 in record_features:
                                                if feature2.type=="CDS":
                                                    CDS_start=np.min([int(feature2.location.start),int(feature2.location.end)])
                                                    CDS_end=np.max([int(feature2.location.start),int(feature2.location.end)])
                                                    if GRINS_start>=(CDS_start-length_var) and GRINS_start<(CDS_end+length_var) and GRINS_end>(CDS_start-length_var) and GRINS_end<=(CDS_end+length_var) and (CDS_end-CDS_start)<500000:
                                                        # GRINS in CDS
                                                        in_CDS=1
                                                        n_GRINS_CDS+=1
                                                        curr_feature=feature2
                                                        break


                                            # recording information about GRINS in the accession-specific output files
                                            GRINS_output_file.write("\t".join([record.id,"GRINSdetect","GRINS",str(GRINS_start),str(GRINS_end),".","+",".",""])+"\n")

                                            if in_BGC==1:
                                                GRINS_BGC_output_file.write("\t".join([record.id,"GRINSdetect","GRINS",str(GRINS_start),str(GRINS_end),".","+",".",curr_BGC_type])+"\n")

                                            if in_CDS==1:
                                                if 'locus_tag' not in curr_feature.qualifiers:
                                                    locus_tag="N/a"
                                                else:
                                                    locus_tag=curr_feature.qualifiers.get('locus_tag')[0]
                                                if 'gene_functions' not in curr_feature.qualifiers:
                                                    gene_function="N/a"
                                                else:
                                                    gene_function=curr_feature.qualifiers.get('gene_functions')[0]
                                                GRINS_CDS_output_file.write("\t".join([record.id,"GRINSdetect","GRINS",str(GRINS_start),str(GRINS_end),".","+",".","locus_tag="+locus_tag+",gene_function="+gene_function])+"\n")

                                # adding the current record to the output GenBank file
                                SeqIO.write(record_annotated, annotated_file, 'genbank')

                                # If specified, plotting the skew plots with duplicated and GRINS regions highlighted
                                if args.with_plots=="yes" or args.with_plots=="Yes" or args.with_plots=="y" or args.with_plots=="Y":

                                    if not isdir(args.plot_output+"/"+assembly):
                                        os.mkdir(args.plot_output+"/"+assembly)

                                    for i in range(0,len(record.seq),100000):

                                        start=i
                                        end=i+100000

                                        curr_duplicated_region_starts=[]
                                        curr_duplicated_region_ends=[]
                                        for item in dupl_locations:
                                            # print item[0],item[1]
                                            if (item[0]>start and item[0]<end) or (item[1]>start and item[1]<end):
                                                # print "bingo!"
                                                curr_duplicated_region_starts.append(np.max([start,item[0]]))
                                                curr_duplicated_region_ends.append(np.min([end,item[1]]))
                                        curr_GRINS_starts=[]
                                        curr_GRINS_ends=[]
                                        for item in GRINS_locations:
                                            if (item[0]>start and item[0]<end) or (item[1]>start and item[1]<end):
                                                curr_GRINS_starts.append(np.max([start,item[0]]))
                                                curr_GRINS_ends.append(np.min([end,item[1]]))

                                        # if len(curr_GRINS_starts)>0:
                                        GCskew_array=GC_skew(str(record.seq[start:end]))
                                        TAskew_array=TA_skew(str(record.seq[start:end]))

                                        plot_graph(args.plot_output,assembly,record.id,start,end,GCskew_array,TAskew_array,curr_duplicated_region_starts,curr_duplicated_region_ends,curr_GRINS_starts,curr_GRINS_ends)

            # Recording information about this accession in the overall output file
            output_file.write("\t".join([assembly,str(length_genome),str(n_contigs_genome),str(n_dups),str(n_BGC),str(length_BGC),str(length_PKS),str(length_NRPS),str(length_terpene),str(length_lanthi),str(length_ladderane),str(length_nucleoside),str(n_GRINS_total),str(n_GRINS_CDS),str(n_GRINS_BGC),str(n_GRINS_PKS),str(n_GRINS_NRPS),str(n_GRINS_terpenes),str(n_GRINS_lanthi),str(n_GRINS_ladderane),str(n_GRINS_nucleoside)])+"\n")
