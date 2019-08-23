#!/usr/bin/env python
# Copyright (C) 2019 Sur Herrera Paredes

# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.

# You should have received a copy of the GNU General Public License
# along with this program.  If not, see <http://www.gnu.org/licenses/>.


import argparse
from Bio import SeqIO
from Bio import SeqRecord
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes a GFF3 file and a sequence file and "
                          "produces a fasta file with one record per "
                          "feature in the GFF3 file.")

    # Define required arguments
    required.add_argument("--input", help=("Input sequence file. Either "
                                           "in fasta on genbank format."),
                          required=True, type=str)
    required.add_argument("--gff3", help=("GFF3 file with header"),
                          required=True,
                          type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Name of the output fasta file. "
                                          "If '', then the input filenames "
                                          "without extension is used and the "
                                          "'.gff3.fasta' suffix is added."),
                        type=str,
                        default='')
    parser.add_argument("--format", help=("Format of the input file."),
                        type=str,
                        default='fasta',
                        choices=['fasta', 'genbank'])

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.output == '':
        args.output = os.path.basename('/home/sur/test.txt')
        args.output = os.path.splitext(args.output)[0]
        args.output = args.output + '.gff3.fasta'
        if os.path.isfile(args.output):
            raise ValueError("{} file already exists.".format(args.output))

    return args


if __name__ == "__main__":
    args = process_arguments()

    seq = SeqIO.read(args.input, args.format)

    # Read windows & write file with pGRINS
    with open(args.gff3, 'r') as ih, open(args.output, 'w') as oh:
        i = 0
        for line in ih:
            if i == 0:
                # Skipping header
                i = i + 1
                continue
            line = line[:-1]
            line = line.split("\t")

            start = int(line[3]) - 1
            end = int(line[4]) - 1
            seqid = line[8].split(';')[0].split('=')[1]
            grins_seq = SeqRecord.SeqRecord(id=seqid,
                                            seq=seq.seq[start:end+1],
                                            description='')
            SeqIO.write(grins_seq, oh, 'fasta')
    ih.close()
    oh.close()
