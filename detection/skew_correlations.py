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
import pybedtools as bed
from Bio import SeqIO
import pandas as pd
import scipy.stats as stats


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
                                          "'_cors.bed'."),
                        type=str,
                        default='')
    # parser.add_argument("--format", help=("format of the input file"),
    #                     default="gff3",
    #                     type=str,
    #                     choices=['gff3'])
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
        args.output = args.output + '_cors.bed'
        if os.path.isfile(args.output):
            raise ValueError("{} file already exists.".format(args.output))

    return args


def interval_correlations(intervals, Skews):
    """Calculate correlations per interval"""

    Res = []
    i = 1
    for w in intervals:
        # print(w)
        ii = (Skews.pos >= w.start) & (Skews.pos <= w.end)
        cor = stats.pearsonr(Skews.GC[ii], Skews.AT[ii])
        if len(w.attrs) == 0:
            ID = "ND_" + str(i)
        else:
            ID = w.attrs['ID']
        Res.append([w.chrom, w.start, w.end, ID, cor[0]])
        i = i + 1

    return(Res)


if __name__ == "__main__":
    args = process_arguments()

    pgrins = bed.BedTool(args.input) \
        .filter(lambda x: len(x) >= args.min_size).sort()
    Skews = pd.read_csv(args.skews,
                        sep="\t",
                        header=None,
                        names=['ref', 'pos', 'GC', 'AT'])
    Res = interval_correlations(intervals=pgrins, Skews=Skews)

    if args.include_complement:
        seqlen = len(SeqIO.read(args.sequence, 'fasta'))
        pgrins.set_chromsizes({'seq': [0, seqlen]})
        non_grins = pgrins.complement() \
            .filter(lambda x: len(x) >= args.min_size).sort()
        Res.extend(interval_correlations(intervals=non_grins, Skews=Skews))

    with open(args.output, 'w') as oh:
        for line in Res:
            line = "\t".join([str(i) for i in line])
            oh.write(line + "\n")
    oh.close()
