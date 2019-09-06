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


from Bio import SeqIO
# import matplotlib.pyplot as plt
# import matplotlib.cm as cm
# import matplotlib.patches as patches
# import numpy as np
# import produce_windows_from_bam as wfb
import argparse
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Calculate sequence skews over sliding windows")

    # Define required arguments
    required.add_argument("--input", help=("Input sequence file."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--prefix", help=("Outfile name prefix. If '' is "
                                          "passed then the name of input "
                                          "without the extension is used."),
                        type=str,
                        default='')
    parser.add_argument("--format", help=("format of the input file"),
                        default="fasta",
                        type=str,
                        choices=['fasta', 'genbank'])
    parser.add_argument("--w_size", help=("Window size."),
                        type=int,
                        default=150)
    parser.add_argument("--s_size", help=("Step size for the sliding window."),
                        type=int,
                        default=30)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.prefix == '':
        args.prefix = os.path.basename(args.input)
        args.prefix = os.path.splitext(args.prefix)[0]

    return args


def gc_skew(seq, w_size=150, s_size=30):
    """Calculate GC Skew"""

    r_len = len(seq)
    Skews = []
    for pos in range(0, r_len - w_size + 1, s_size):
        end = min(pos+w_size, r_len)
        window = seq[pos:end]
        gc_skew = ((window.seq.count('G') - window.seq.count('C'))
                   / (window.seq.count('G') + window.seq.count('C')))
        Skews.append([pos, end, gc_skew])

    return Skews


def at_skew(seq, w_size=150, s_size=30):
    """Calculate AT Skew"""

    r_len = len(seq)
    Skews = []
    for pos in range(0, r_len - w_size + 1, s_size):
        end = min(pos + w_size, r_len)
        window = seq[pos:end]
        at_skew = ((window.seq.count('A') - window.seq.count('T'))
                   / (window.seq.count('A') + window.seq.count('T')))
        Skews.append([pos, end, at_skew])

    return Skews


if __name__ == "__main__":
    args = process_arguments()

    seq = SeqIO.read(args.input, args.format)
    gc_skew = gc_skew(seq, w_size=args.w_size, s_size=args.s_size)
    at_skew = at_skew(seq, w_size=args.w_size, s_size=args.s_size)

    # Write output, assume identical lengths of skew lists
    outfile = args.prefix + "_skews.txt"
    with open(outfile, 'w') as oh:
        for i in range(len(gc_skew)):
            pos_gc = (gc_skew[i][0] + gc_skew[i][1]) / 2
            pos_at = (at_skew[i][0] + at_skew[i][1]) / 2
            if pos_at != pos_gc:
                raise ValueError("Skew positions don't match")
            line = [seq.id, str(pos_gc),
                    str(gc_skew[i][2]),
                    str(at_skew[i][2])]
            line = "\t".join(line)
            oh.write(line + "\n")
