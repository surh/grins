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
import matplotlib.pyplot as plt
import matplotlib.cm as cm
import matplotlib.patches as patches
import numpy as np
import produce_windows_from_bam as wfb
import argparse
import os


def gc_skew(seq, w_size=150, s_size=30):
    """Calculate GC Skew"""

    r_len = len(seq)
    Skews = []
    for pos in range(0, r_len - w_size + 1, s_size):
        end = min(pos+w_size, r_len)
        window = seq[pos:end]
        denom = window.seq.count('G') + window.seq.count('C')
        if denom == 0:
            denom = 1
        gc_skew = ((window.seq.count('G') - window.seq.count('C'))
                   / denom)
        Skews.append([pos, end, gc_skew])

    return Skews


def read_gff3(gff3_file):
    """Read windows from a GFF3 file with header."""

    # with open(windows_file, 'r') as ih, open(pgrins_file, 'w') as oh:
    with open(gff3_file, 'r') as ih:
        i = 0
        p_grins = []
        for line in ih:
            if i == 0:
                i = i + 1
                continue
            line = line[:-1]
            line = line.split("\t")
            p_grins.append(line)
    ih.close()

    return p_grins


def read_uc(uc_file):
    """Read cluster assignments from .uc file"""

    with open(uc_file, 'r') as ih:
        Clusters = dict()
        n_clusters = 0
        for line in ih:
            line = line.split("\t")
            Clusters[line[8]] = int(line[1])
            if int(line[1]) > n_clusters:
                n_clusters = int(line[1])
        n_clusters = n_clusters + 1
    ih.close()

    return Clusters, n_clusters


def plot_fast_grins(p_grins, Skews, bam_windows,
                    Clusters, outfile, n_clusters,
                    w_size=150):
    """Plot results of fast grins prediction"""

    plt.figure(figsize=(20, 10), dpi=150)
    y_max = np.max([i[2] for i in Skews])
    y_min = np.min([i[2] for i in Skews])
    new_max = y_max + (y_max - y_min) / 4
    plt.ylim(top=new_max, bottom=y_min)
    for g in p_grins:
        start = int(g[3])
        end = int(g[4])
        grins_id = g[8].split(';')[0].split('=')[1]
        breaks = np.linspace(0, 1, n_clusters)
        plt.axvspan(start, end,
                    facecolor=cm.gist_rainbow(breaks)[Clusters[grins_id]],
                    alpha=0.8)
    midpoints = [(x[0] + x[1]) / 2 for x in Skews]
    skews = [x[2] for x in Skews]
    plt.plot(midpoints, skews, color='black')

    conn_style = patches.ConnectionStyle.Arc3(rad=-0.5)
    for w in bam_windows:
        w_pos = w[1] + (w_size / 2)
        m_pos = w[3] + (w_size / 2)
        pos = np.sort([w_pos, m_pos])
        p = patches.FancyArrowPatch(posA=(pos[0], y_max),
                                    posB=(pos[1], y_max),
                                    connectionstyle=conn_style,
                                    arrowstyle='-',
                                    color='grey',
                                    linewidth=1)
        plt.gca().add_patch(p)

    plt.savefig(outfile)


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Plot the results of detect_fast for grins")

    # Define required arguments
    required.add_argument("--input", help=("Input sequence file."),
                          required=True, type=str)
    required.add_argument("--grins_gff3",
                          help=("A GFF3 file containing only GRINS."),
                          required=True,
                          type=str)
    required.add_argument("--grins_clusters",
                          help=("A .uc file produced by vsearch, "
                                "clustering the GRINS in --grins_gff3."),
                          required=True,
                          type=str)
    required.add_argument("--windows_bam",
                          help=("A BAM file produced by bowtie2 with the "
                                "result of mapping every window in a sequence "
                                "to the rest of the sequence"),
                          required=True,
                          type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Outfile name. If '' is passed "
                                          "then the name of input with the "
                                          "extension changed for .png."),
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
    if args.output == '':
        args.output = os.path.basename(args.input)
        args.output = os.path.splitext(args.output)[0]
        args.output = args.output + '.png'
        if os.path.isfile(args.output):
            raise ValueError("{} file already exists.".format(args.output))

    return args


if __name__ == "__main__":
    args = process_arguments()

    seq = SeqIO.read(args.input, args.format)
    Skews = gc_skew(seq, w_size=args.w_size, s_size=args.s_size)
    p_grins = read_gff3(args.grins_gff3)
    Clusters, n_clusters = read_uc(args.grins_clusters)
    bam_windows = wfb.find_bam_windows(file=args.windows_bam)
    plot_fast_grins(p_grins=p_grins,
                    Skews=Skews,
                    bam_windows=bam_windows,
                    Clusters=Clusters,
                    outfile=args.output,
                    n_clusters=n_clusters,
                    w_size=args.w_size)
