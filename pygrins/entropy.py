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


def calculate_compression_ratio(seq):
    seq_c = zlib.compress(seq.encode('utf-8'), level=9)

    return (len(seq_c) - 8) / len(seq)


def sw_compression_ratio(seq, start=0, step=1, window=1):
    H_seq = []
    for w_start in range(start, len(seq) - window + 1, step):
        w_end = min(w_start + window, len(seq))
        H_seq.append([w_start, w_end,
                      calculate_compression_ratio(seq[w_start:w_end])])

    return H_seq


def kmer_shannon(seq, k=3):
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
    """Sliding window k-mer based Shannon entropy. Avoids re-calculating
    shared sections between windows"""

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
