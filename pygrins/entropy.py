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

import zlib
import numpy as np
import argparse


def calculate_compression_ratio(seq):
    """Compresses a string using zlib and returns the compression
    ratio, as the ratio of number of characters in the compressed
    string over the number of characters in the original screen

    :param seq: A sequence represented as a string.
    :type seq: str

    :return: The compression ratio of the sequence by zlib
    :rtype: float"""

    if type(seq) is not str:
        raise ValueError("seq must be of type 'str'.")

    seq_c = zlib.compress(seq.encode('utf-8'), level=9)

    # The minus 8 is the decompression key string length.
    return (len(seq_c) - 8) / len(seq)


def sw_compression_ratio(seq, start=0, step=1, window=1):
    """Takes a sequence, splits it in a series of sliding windows,
    and applies calculate_compression_ratio to each window.

    :param seq: A sequence represented as a string
    :type seq: str
    :param start: Position to start the slidinw windows (0-indexed).
    :type start: int, optional
    :param step: Step size of the sliding window.
    :type step: int, optional
    :param window: Window size of each sliding window.
    :type window: int, optional

    :return: A list of lists where each element gives the start and
    end positions (0-indexed, open interval on the right) of each window,
    and the compression ratio of each window.
    :rtype: list
    """

    if start < 0:
        raise ValueError("start must be a non-negative integer.")
    if step < 1 or window < 1:
        raise ValueError("step and window must be positive integers.")

    H_seq = []
    for w_start in range(start, len(seq) - window + 1, step):
        w_end = min(w_start + window, len(seq))
        H_seq.append([w_start, w_end,
                      calculate_compression_ratio(seq[w_start:w_end])])

    return H_seq


def kmer_shannon(seq, k=3):
    """Calculates the k-mer based Shannon entropy of a given sequence,
    and k-mer length k.

    :param seq: A sequence represented as a string.
    :type seq: str
    :param k: Length for the k-mers.
    :type k: int, optional

    :return: The k-mer based Shannon entropy of the sequence.
    :rtype: float"""

    if type(seq) is not str:
        raise ValueError("seq must be of type 'str'.")
    if k < 1:
        raise ValueError("k must be a be positive integer.")

    Kmers = dict()
    for start in range(len(seq) - k + 1):
        kmer = seq[start:start+k]
        if kmer in Kmers:
            Kmers[kmer] = Kmers[kmer] + 1
        else:
            Kmers[kmer] = 1

    counts = np.fromiter(Kmers.values(), dtype=int)
    pi = counts / np.sum(counts)
    Hk_seq = -np.sum(pi * np.log2(pi))

    return Hk_seq


def sw_kmer_shannon(seq, start=0, step=5, window=20, k=3):
    """ Takes a sequence, splits it in a series of sliding windows,
    and calculates the k-mer based Shannon entropy. It is more
    efficient than independently applying kmer_shannon to each window,
    because it avoids re-calculating the kmer abundances for shared portions
    of consecutive sliding windows.

    :param seq: A sequence represented as a string.
    :type seq: str
    :param start: Position to start the slidinw windows (0-indexed).
    :type start: int, optional
    :param step: Step size of the sliding window.
    :type step: int, optional
    :param window: Window size of each sliding window.
    :type window: int, optional
    :param k: Length for the k-mers.
    :type k: int, optional

    :return: A list of lists where each element gives the start and
    end positions (0-indexed, open interval on the right) of each window,
    and the compression ratio of each window.
    :rtype: list
    """

    if type(seq) is not str:
        raise ValueError("seq must be of type 'str'.")
    if step < 1 or window < 1 or k < 1:
        raise ValueError("step, window and k must be positive integers.")

    prev_kmers = []
    Hk_seq = []
    Kmers = dict()
    initialized = False
    for w_start in range(start, len(seq) - window + 1, step):
        w_end = min(w_start + window, len(seq))
        fast_start = w_start

        if initialized:
            fast_start = np.maximum(w_end - step - k + 1, w_start)
            for kmer in prev_kmers[0:step]:
                Kmers[kmer] = Kmers[kmer] - 1
                # print("\t=rm", kmer)
            prev_kmers = prev_kmers[step:]

        # print(w_start, fast_start, w_end)
        for start in range(fast_start, w_end - k + 1):
            kmer = seq[start:start+k]
            # print("\t", start, start+k, kmer)

            prev_kmers.append(kmer)

            if kmer in Kmers:
                Kmers[kmer] = Kmers[kmer] + 1
            else:
                Kmers[kmer] = 1
                initialized = True
        w_counts = np.fromiter(Kmers.values(), dtype=int)
        w_counts = w_counts[w_counts != 0]
        w_pi = w_counts / np.sum(w_counts)
        Hk_seq.append([w_start, w_end, -np.sum(w_pi * np.log2(w_pi))])

    return Hk_seq


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that takes a fasta file and, for each "
                          "record in the file, it calculates its compression "
                          "ratio, and k-mer based Shannon entropy for various "
                          "k values. It prints the results into a file. "
                          "The functions in the file, can also be imported "
                          "and used in another script. It requires numpy and "
                          "zlib. All cordinates are 0-indexed.")

    # Define required arguments
    required.add_argument("fasta", help=("FASTA file with sequences."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outfile", help=("Name of the outfile for results"),
                        type=str,
                        default="results.txt")
    parser.add_argument("--step", help=("Step size."),
                        default=30,
                        type=int)
    parser.add_argument("--window", help=("Window size."),
                        type=int,
                        default=150)
    parser.add_argument("-k", help=("Sizes for k-mers to calculate entropy."),
                        type=int,
                        default=[3, 4])
    parser.add_argument("--format", help=("Format of output table."),
                        type=str,
                        default="bed",
                        choices=["bed", "gff", "tab"])

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    print(args.k)
