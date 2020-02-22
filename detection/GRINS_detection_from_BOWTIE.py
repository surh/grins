import Bio
from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.SeqFeature import SeqFeature,FeatureLocation
import numpy as np
import argparse
import os
from os import listdir
from os.path import isfile, isdir, join
# import matplotlib
# matplotlib.use('TkAgg')
# import matplotlib.pyplot as plt
# import matplotlib.patches as patches
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

    required.add_argument("--duplBGC_input", help=("Folder with BGC duplication detection results"),
                          required=True, type=str)

    parser.add_argument("--seq_output", help=("Annotated GenBank output folder"),
                          required=False, type=str, default="./output/genomes_GRINS")

    parser.add_argument("--GRINS_output", help=("GRINS output folder"),
                          required=False, type=str, default="./output/GRINS.gff3")

    parser.add_argument("--GRINS_BGC_output", help=("GRINS in BGC output folder"),
                          required=False, type=str, default="./output/GRINS_BGC.gff3")
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

# # Plotting the graph for each record
# def plot_graph(name,start,end,GCskew_array,TAskew_array,dupl_region_starts,dupl_region_ends,GRINS_starts,GRINS_ends):
#         a=np.arange(start,start+len(GCskew_array)*30,30)
#         plt.figure(0)
#         fig,ax1=plt.subplots(figsize=(10, 4))
#         ax1.set_ylim(-1,1)
#         ax1.set_xlim(start,end)
#         ax1.grid(False)
#         plt.plot(a,GCskew_array,color="blue",linewidth=0.5)
#         plt.plot(a,TAskew_array,color="red",linewidth=0.5)
#         plt.legend(["GC skew","TA skew"],loc=4)
#         plt.legend(["Pairwise similarity"],loc=3)
#         for j in range(0,len(dupl_region_starts)):
#             region_start=dupl_region_starts[j]
#             region_end=dupl_region_ends[j]
#             rect = patches.Rectangle((region_start,1),region_end-region_start,-2,linewidth=1,edgecolor='grey',facecolor='grey', alpha=0.3)
#             ax1.add_patch(rect)
#         for k in range(0,len(GRINS_starts)):
#             region_start=GRINS_starts[k]
#             region_end=GRINS_ends[k]
#             rect = patches.Rectangle((region_start,1),region_end-region_start,-2,linewidth=1,edgecolor='teal',facecolor='teal', alpha=0.3)
#             ax1.add_patch(rect)


#         plt.title("GRINS in %s, %d-%d kb" % (name,start/1000,end/1000))
#         plt.savefig("./%s/Plots/%s_GRINS_%d-%dkb.png" % (name,name,start/1000,end/1000) )
#         plt.close('all')

if __name__ == "__main__":

    args = process_arguments()
    assembly_accessions=[f[:f.find("_genomic.duplicated.gff3")] for f in listdir(args.dupl_input) if isfile(join(args.dupl_input, f))]
    sequence_files=[f for f in listdir(args.seq_input) if isfile(join(args.seq_input, f))]

    if not isdir(args.seq_output):
        os.mkdir(args.seq_output)
    if not isdir(args.GRINS_output):
        os.mkdir(args.GRINS_output)
    if not isdir(args.GRINS_BGC_output):
        os.mkdir(args.GRINS_BGC_output)

    for assembly in assembly_accessions:

        # Reading duplication detection results
        with open("%s/%s_genomic.duplicated.gff3" %(args.dupl_input,assembly),'r') as dupl_results_file:
            dupl_results=dupl_results_file.readlines()
        dupl_location_dict={}
        for i in range(1,len(dupl_results)):
            line=dupl_results[i].strip("\n").split("\t")
            accession=line[0]
            start=int(line[3])
            end=int(line[4])
            if accession in dupl_location_dict:
                dupl_location_dict[accession].append([start,end])
            else:
                dupl_location_dict[accession]=[[start,end]]

        # Reading duplication in BGC detection results
        with open("%s/%s_genomic.bgcdups.gff3" %(args.duplBGC_input,assembly),'r') as dupl_BGC_results_file:
            dupl_BGC_results=dupl_BGC_results_file.readlines()
        dupl_BGC_location_dict={}
        for i in range(0,len(dupl_BGC_results)):
            line=dupl_BGC_results[i].strip("\n").split("\t")
            accession=line[0]
            start=int(line[3])
            end=int(line[4])
            attributes=line[8]
            if accession in dupl_BGC_location_dict:
                dupl_BGC_location_dict[accession].append([start,end,attributes])
            else:
                dupl_BGC_location_dict[accession]=[[start,end,attributes]]

        # Reading sequence file
        filename = [s for s in sequence_files if assembly in s]
        records=SeqIO.parse("./%s/%s" %(args.seq_input,filename[0]),'gb')
        length_var=300

        # Writing a new record, with detected GRINS annotated
        with open("./%s/%s.GRINS.gbk" %(args.seq_output,assembly),"w") as annotated_file:

            with open("./%s/%s.GRINS_BGC.gff3" %(args.GRINS_BGC_output,assembly),'w') as GRINS_BGC_output_file:
                # GRINS_output_file.write("\t".join(["Record","GRINS start","GRINS end","Locus","CDS name"])+"\n")
                
                with open("./%s/%s.GRINS.gff3" %(args.GRINS_output,assembly),'w') as GRINS_output_file:
                    # GRINS_nonbiosynth_output_file.write("\t".join(["Record","GRINS start","GRINS end","Locus","CDS name"])+"\n")

                    for record in records:

                        record_features=record.features
                        record_annotated=record
                        accession=record.id
                        GRINS_locations=[]
                        if accession in dupl_location_dict:
                            dupl_locations=dupl_location_dict[accession]
                        else:
                            dupl_locations=[]

                        if accession in dupl_BGC_location_dict:
                            dupl_BGC_locations=[[i[0],i[1]] for i in dupl_BGC_location_dict[accession]]
                            dupl_BGC_clusters=[i[2] for i in dupl_BGC_location_dict[accession]]
                        else:
                            dupl_BGC_locations=[]
                            dupl_BGC_clusters=[]

                        for i in range(0,len(dupl_locations)):
                            if dupl_locations[i][1]-dupl_locations[i][0]>=500:
                                mean_GCskew,mean_TAskew=abs_skew_means(record,dupl_locations[i][0],dupl_locations[i][1],150,30)
                                if mean_GCskew>=0.15 and mean_TAskew>=0.15:
                                    in_CDS=0
                                    in_biosynth_CDS=0
                                    feature=SeqFeature(FeatureLocation(start=dupl_locations[i][0], end=dupl_locations[i][1]), type='GRINS')
                                    record_annotated.features.append(feature)
                                    GRINS_start=int(dupl_locations[i][0])
                                    GRINS_end=int(dupl_locations[i][1])
                                    GRINS_locations.append((GRINS_start,GRINS_end))

                                    # checking if GRINS in a CDS
                                    for feature1 in record_features:
                                        if feature1.type=="CDS":
                                            CDS_start=np.min([int(feature1.location.start),int(feature1.location.end)])
                                            CDS_end=np.max([int(feature1.location.start),int(feature1.location.end)])
                                            if GRINS_start>=(CDS_start-length_var) and GRINS_start<(CDS_end+length_var) and GRINS_end>(CDS_start-length_var) and GRINS_end<=(CDS_end+length_var) and (CDS_end-CDS_start)<500000:
                                                in_CDS=1
                                                curr_feature=feature1

                                                # checking if GRINS in BGC
                                                if dupl_locations[i] in dupl_BGC_locations:
                                                    index=dupl_BGC_locations.index(dupl_locations[i])
                                                    curr_clutster=dupl_BGC_clusters[index]
                                                    in_biosynth_CDS=1

                                                break

                                    if in_biosynth_CDS==1:
                                        # if curr_feature.qualifiers.get('product')[0]
                                        if 'locus_tag' not in curr_feature.qualifiers:
                                            locus_tag="N/a"
                                        else:
                                            locus_tag=curr_feature.qualifiers.get('locus_tag')[0]
                                        if 'gene_functions' not in curr_feature.qualifiers:
                                            gene_function="N/a"
                                        else:
                                            gene_function=curr_feature.qualifiers.get('gene_functions')[0]
                                        GRINS_output_file.write("\t".join([record.id,"GRINSdetect","GRINS",str(GRINS_start),str(GRINS_end),".","+",".","locus_tag="+locus_tag+",gene_function="+gene_function+","+curr_clutster])+"\n")
                                        GRINS_BGC_output_file.write("\t".join([record.id,"GRINSdetect","GRINS",str(GRINS_start),str(GRINS_end),".","+",".","locus_tag="+locus_tag+",gene_function="+gene_function+","+curr_clutster])+"\n")

                                    elif in_CDS==1:
                                        if 'locus_tag' not in curr_feature.qualifiers:
                                            locus_tag="N/a"
                                        else:
                                            locus_tag=curr_feature.qualifiers.get('locus_tag')[0]
                                        if 'gene_functions' not in curr_feature.qualifiers:
                                            gene_function="N/a"
                                        else:
                                            gene_function=curr_feature.qualifiers.get('gene_functions')[0]
                                        GRINS_output_file.write("\t".join([record.id,"GRINSdetect","GRINS",str(GRINS_start),str(GRINS_end),".","+",".","locus_tag="+locus_tag+",gene_function="+gene_function])+"\n")
                                    else:
                                        GRINS_output_file.write("\t".join([record.id,"GRINSdetect","GRINS",str(GRINS_start),str(GRINS_end),".","+",".","locus_tag=outside_CDS"])+"\n")


                            feature=SeqFeature(FeatureLocation(start=dupl_locations[i][0], end=dupl_locations[i][1]), type='Duplication')
                            record_annotated.features.append(feature)            


                        SeqIO.write(record_annotated, annotated_file, 'genbank')

                        # for i in range(0,len(sequence),100000):
                        # # for i in range(1120000,1220000,100000):
                        #     start=i
                        #     end=i+100000
                        #     # print ""
                        #     # print start,end
                        #     GCskew_array=GC_skew(str(sequence[start:end]))
                        #     TAskew_array=TA_skew(str(sequence[start:end]))

                        #     curr_duplicated_region_starts=[]
                        #     curr_duplicated_region_ends=[]
                        #     for item in dupl_locations:
                        #         # print item[0],item[1]
                        #         if (item[0]>start and item[0]<end) or (item[1]>start and item[1]<end):
                        #             # print "bingo!"
                        #             curr_duplicated_region_starts.append(np.max([start,item[0]]))
                        #             curr_duplicated_region_ends.append(np.min([end,item[1]]))
                        #     curr_GRINS_starts=[]
                        #     curr_GRINS_ends=[]
                        #     for item in GRINS_locations:
                        #         if (item[0]>start and item[0]<end) or (item[1]>start and item[1]<end):
                        #             curr_GRINS_starts.append(np.max([start,item[0]]))
                        #             curr_GRINS_ends.append(np.min([end,item[1]]))

                            
                        #     # print "duplications",curr_duplicated_region_starts,curr_duplicated_region_ends
                        #     # print "GRINS",curr_GRINS_starts,curr_GRINS_ends
                        #     plot_graph(species_name,start,end,GCskew_array,TAskew_array,curr_duplicated_region_starts,curr_duplicated_region_ends,curr_GRINS_starts,curr_GRINS_ends)
