#!/usr/bin/env python
# Copyright (C) 2019 Aleksandra Nivina, Sur Herrera Paredes

# import Bio
# from Bio.Seq import Seq
from Bio import SeqIO
# from Bio.Alphabet import IUPAC
from Bio import pairwise2
# import numpy as np
import argparse


def process_arguments():
    """ Read and process"""

    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("User-specified parameters for GRINS detection")

    # Define required arguments
    required.add_argument('-infile ', "--input_file", type=str,
                          help="name and location of the input GenBank file",
                          required=True)

    # Define other arguments
    parser.add_argument('-outfolder ', "--output_folder", type=str,
                        help="location of the output folder",
                        default='./')

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Read sequence. NOTE: it only handles one record per file.
    record = SeqIO.read(args.input_file, "gb")
    sequence = record.seq
    record_name = record.name

    # performing pairwise alignments for a sliding window of
    # 150nt and a step of 30nt
    print("Starting homology matrix")
    homology_matrix = []

    for i in range(0, len(sequence)):
        homology_matrix.append([])
        print("Current window:", i)
        for j in range(0, len(sequence)):
            subseq1 = sequence[i:i + 150]
            subseq2 = sequence[j:j + 150]
            alignment = pairwise2.align.localxx(subseq1, subseq2,
                                                score_only=True)
            homology_matrix[len(homology_matrix)-1].append(alignment)

    # saving the results
    with open(args.output_folder+record_name+"_window150_step30.txt", 'w') as output_file:
        for i in range(0, len(homology_matrix)):
            for j in range(0, len(homology_matrix[i])):
                output_file.write(str(homology_matrix[i][j])+",")
            output_file.write("\n")
