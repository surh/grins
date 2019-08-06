#!/usr/bin/env python
# Copyright (C) 2019 Aleksandra Nivina, Sur Herrera Paredes

from Bio import SeqIO
from Bio import pairwise2
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
    required.add_argument("--input", type=str,
                          help="name and location of the input GenBank file",
                          required=True)

    # Define other arguments
    parser.add_argument("--output", type=str,
                        help="location of the output folder",
                        default='./')
    parser.add_argument("--w_size", help=("Window size in nucleotides"),
                        type=int,
                        default=150)
    parser.add_argument("--s_size", help=("Step size in nucleotides."),
                        type=int,
                        default=30)
    parser.add_argument()

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Read sequence. NOTE: it only handles one record per file.
    record = SeqIO.read(args.input, "gb")
    sequence = record.seq
    record_name = record.name

    # performing pairwise alignments for a sliding window
    print("Initializing homology matrix")
    homology_matrix = []

    for i in range(0, len(sequence), args.s_size):
        homology_matrix.append([])
        print("Current window:", i)
        subseq1 = sequence[i:i + args.w_size]
        for j in range(0, len(sequence), args.s_size):
            subseq2 = sequence[j:j + args.w_size]
            alignment = pairwise2.align.localxx(subseq1, subseq2,
                                                score_only=True)
            homology_matrix[len(homology_matrix)-1].append(alignment)

    # saving the results
    output_name = "_window" + str(args.w_size) \
                  + "_step" + str(args.s_size) + ".txt"
    output_name = args.output + '/' + record_name + output_name
    print("Writing output file")
    with open(output_name, 'w') as output_file:
        for i in range(0, len(homology_matrix)):
            for j in range(0, len(homology_matrix[i])):
                output_file.write(str(homology_matrix[i][j])+",")
            output_file.write("\n")
