import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import pairwise2
import numpy as np
import argparse


parser=argparse.ArgumentParser(description="User-specified parameters for GRINS detection")
# Input file:
parser.add_argument('-infile ',"--input_file",type=str, help="name and location of the input GenBank file")
# Output folder:
parser.add_argument('-outfolder ',"--output_folder",type=str, help="location of the output folder")


args=parser.parse_args()

record = SeqIO.read(args.input_file, "gb")
sequence=record.seq
record_name=record.name

# performing pairwise alignments for a sliding window of 150nt and a step of 30nt
homology_matrix=[]
for i in range(0,len(sequence)):
	homology_matrix.append([])
	for j in range(0,len(sequence)):
		subseq1=sequence[i:i+150]
		subseq2=sequence[j:j+150]
		alignment=pairwise2.align.localxx(subseq1, subseq2, score_only=True)
		homology_matrix[len(homology_matrix)-1].append(alignment)

# saving the results
with open(args.output_folder+record_name+"_window150_step30.txt",'w') as output_file:
    for i in range(0,len(homology_matrix)):
        for j in range(0,len(homology_matrix[i])):
            output_file.write(str(homology_matrix[i][j])+",")
        output_file.write("\n")
     