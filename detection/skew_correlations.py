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

import os
import argparse
import plot_fast_grins as pfg
import pybedtools as bed
from Bio import SeqIO


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Read a skew files and a file with intervals "
                          "and calculate correlations per interval,")

    # Define required arguments
    required.add_argument("--input", help=("File with genomic intervals"),
                          required=True, type=str)
    required.add_argument("--skews", help=("File with skews"),
                          required=True, type=str)
    required.add_argument("--sequence", help=("Fasta file with sequence"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Outfile name. If '' is passed "
                                          "then the name of input with the "
                                          "extension changed for "
                                          "'_cors.txt'."),
                        type=str,
                        default='')
    parser.add_argument("--format", help=("format of the input file"),
                        default="gff3",
                        type=str,
                        choices=['gff3'])
    parser.add_argument("--include_complement",
                        help=("If passed the complement of the intervals "
                              "is calculated and correlations for the "
                              "complement is also returned"),
                        action="store_true",
                        default=False)
    parser.add_argument("--min_size",
                        help=("Minimum interval size"),
                        type=int,
                        default=500)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.output == '':
        args.output = os.path.basename(args.input)
        args.output = os.path.splitext(args.output)[0]
        args.output = args.output + '_cors.txt'
        if os.path.isfile(args.output):
            raise ValueError("{} file already exists.".format(args.output))

    return args


def gff3_to_bed(gff3, outfile=False):
    """Convert a GFF3 list into a BED file"""


if __name__ == "__main__":
    args = process_arguments()

    if args.format == 'gff3':
        gff3 = pfg.read_gff3(args.input)
