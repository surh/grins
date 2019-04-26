from Bio import SeqIO
from Bio.Seq import Seq
from Bio.SeqRecord import SeqRecord
from Bio.Alphabet import IUPAC
import os
from os import listdir
from os.path import isfile, join
import argparse


def cluster_number(cluster_n):
    if len(cluster_n)==1:
        cluster_num="00"+cluster_n
    elif len(cluster_n)==2:
        cluster_num="0"+cluster_n
    else:
        cluster_num=cluster_n
    return cluster_num


parser=argparse.ArgumentParser(description="Location of the input files")

parser.add_argument('-infolder ',"--input_folder",type=str, help="Location of the folder with input GenBank files")
parser.add_argument('-known ',"--known_cluster_file",type=str, help="Location of the file with a list of known clusters")
parser.add_argument('-outfolder ',"--output_folder",type=str, help="Location of the output folder")
args=parser.parse_args()

# Reading the list of known clusters
with open(args.known_cluster_file,'r') as known_file:
    known_cluster_list=known_file.readlines()


# For each known cluster, extracting annotated biosynthetic CDSs and saving tehm into a GenBank file
for known_cluster in known_cluster_list:
    known_cluster_data=known_cluster.strip("\n").split("\t")
    known_cluster_id=known_cluster_data[0].strip(":").split(" ")
    accession=known_cluster_id[0]
    cluster_n=known_cluster_id[1]

    if accession=="AY77199" and cluster_n=="1":
        accession="AM420293"
        cluster_n=="2"

    cluster_num=cluster_number(cluster_n)
    my_in_path=args.input_folder+accession+"/"
    file=accession+".cluster"+cluster_num+".gbk"

    print ("Working on",my_in_path+file)
    record = SeqIO.read(my_in_path+file, "gb")
    feature_list=record.features
    feature_n=len(feature_list)
    new_record=SeqRecord("",id=record.id,name=record.name,description=record.description)

    with open(args.output_folder+file, "w") as output_file:   

        # Counting the strand with more ORFs
        forward_strand_CDS=0
        reverse_strand_CDS=0
        CDS_orientation=""

        # going over CDS features
        for n1 in range(0,feature_n):
            feature1=feature_list[n1]
            if feature1.type=="CDS":

                # going over aSDomain features within that CDS
                for n2 in range(n1+1,feature_n):
                    feature2=feature_list[n2]
                    if feature2.location.start in feature1 and (feature2.type=="aSDomain"):
                        if feature1.location.strand==1:
                            forward_strand_CDS+=1
                            CDS_orientation+="F"
                        else:
                            reverse_strand_CDS+=1
                            CDS_orientation+="R"
                        break


        # If more ORFs in the Forward direction
        if forward_strand_CDS>=reverse_strand_CDS:
            # going over CDS features
            current_record=SeqRecord("",id=record.id,name=record.name,description=record.description)
            prev_orientation=1
            for i1 in range(0,feature_n):
                feature1=feature_list[i1]
                if feature1.type=="CDS":

                    # going over aSDomain features within that CDS
                    for i2 in range(i1+1,feature_n):
                        feature2=feature_list[i2]
                        if feature2.location.start in feature1 and (feature2.type=="aSDomain"):

                            # If located on the Fw strand:
                            if feature1.location.strand==1:
                                additional_record=record[feature1.location.start:feature1.location.end]

                                # Record current (=previous) ORFs
                                new_record=new_record+current_record

                                # Set this ORF as current
                                current_record=additional_record
                                prev_orientation=1
                                break
                            
                            # If located on the Rev strand:
                            else:
                                additional_record=record[feature1.location.start:feature1.location.end].reverse_complement()
                                
                                # If previous was located on the Fw strand:
                                if prev_orientation==1:

                                    # Record current (=previous) ORFs
                                    new_record=new_record+current_record

                                    # Set this ORF as current
                                    current_record=additional_record
                                    prev_orientation=-1
                                    break
                                
                                # If both previous and current ORFs are on the Rev strand, 
                                else:

                                    # Put this and the current (=previous) ORF(s) in opposite order but not record them for now
                                    current_record=additional_record+current_record
                                    prev_orientation=-1
                                    break

            # recording last ORF/set of ORFs
            new_record=new_record+current_record


        else:
            # more ORFs in the Reverse direction
            current_record=SeqRecord("",id=record.id,name=record.name,description=record.description)
            prev_orientation=1
            for i1 in range(feature_n-1,-1,-1):
                feature1=feature_list[i1]
                if feature1.type=="CDS":

                    # going over aSDomain features within that CDS
                    for i2 in range(i1+1,feature_n):
                        feature2=feature_list[i2]
                        if feature2.location.start in feature1 and (feature2.type=="aSDomain"):

                            # If located on the Rev strand:
                            if feature1.location.strand==-1:

                                additional_record=record[feature1.location.start:feature1.location.end].reverse_complement()

                                # Record current (=previous) ORF
                                new_record=new_record+current_record

                                # Set this ORF as current
                                current_record=additional_record
                                prev_orientation=-1
                                break
                            
                            # If located on the Fw strand:
                            else:

                                additional_record=record[feature1.location.start:feature1.location.end]

                                # If previous was located on the Rev strand:
                                if prev_orientation==-1:

                                    # Record current (=previous) ORFs
                                    new_record=new_record+current_record

                                    # Set this ORF as current
                                    current_record=additional_record
                                    prev_orientation=1
                                    break
                                
                                # If previous was located on the Fw strand:
                                else:

                                    # Put this and the current (=previous) ORF(s) in opposite order but not record them for now
                                    current_record=additional_record+current_record
                                    prev_orientation=1
                                    break
            
            # recording last ORF/set of ORFs
            new_record=new_record+current_record
        
        # save this record into a file
        new_record.id=record.id
        new_record.name=record.name
        new_record.description=record.description
        SeqIO.write(new_record, output_file, "gb")


        
        