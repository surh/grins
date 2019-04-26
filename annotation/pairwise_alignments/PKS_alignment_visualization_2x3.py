import Bio
from Bio.Seq import Seq
from Bio import SeqIO
from Bio.Alphabet import IUPAC
from Bio import pairwise2
import numpy as np
import argparse
import time

print "alignment 2x3 started"

parser=argparse.ArgumentParser(description="User-specified parameters for GRINS detection")
# Input file:
parser.add_argument('-infile1 ',"--input_file1",type=str, help="name and location of the input GenBank file")
# Second input file (optional):
parser.add_argument('-infile2 ',"--input_file2",type=str, help="name and location of the input GenBank file")
# Output folder:
parser.add_argument('-outfolder ',"--output_folder",type=str, help="location of the output folder")
# Size of split fragments:
parser.add_argument('-split1 ',"--splitsize1",type=int, help="size of each split gene size (typically ~1/3 of totoal)")
# Size of split fragments for second cluster (optional):
parser.add_argument('-split2 ',"--splitsize2",type=int, help="size of each split gene size (typically ~1/3 of totoal)")

args=parser.parse_args()

record = SeqIO.read(args.input_file1, "gb")
sequence1=record.seq
record_name1=record.name
if args.input_file2!=None:
    record2 = SeqIO.read(args.input_file2, "gb")
    sequence2=record2.seq
    record_name2=record2.name
else:
    sequence2=record.seq
    record_name2=record.name

if args.splitsize2!=None:
    splitsize2=args.splitsize2
else:
    splitsize2=args.splitsize1


homology_matrix=[]
for i in range(args.splitsize1,args.splitsize1*2,30):
	homology_matrix.append([])
	for j in range(splitsize2*2,len(sequence2),30):
		subseq1=sequence1[i:i+150]
		subseq2=sequence2[j:j+150]
		alignment=pairwise2.align.localxx(subseq1, subseq2, score_only=True)
		homology_matrix[len(homology_matrix)-1].append(alignment)

homology_matrix2=np.array(homology_matrix)


all_min=50
all_max=0
for j in range(0,len(homology_matrix)):
	new_min=min(homology_matrix[j])
	if new_min<all_min:
		all_min=new_min
	new_max=max(homology_matrix[j])
	if new_max>all_max:
		all_max=new_max
        
print "min", all_min
print "max", all_max

with open(args.output_folder+record_name1+"_vs_"+record_name2+"_window150_step30_2x3.txt",'w') as output_file:
    for i in range(0,len(homology_matrix)):
        for j in range(0,len(homology_matrix[i])):
            output_file.write(str(homology_matrix[i][j])+",")
        output_file.write("\n")

print "alignment 2x3 finished" 
print time.ctime()       