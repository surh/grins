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




def sw_kmer_shannon(seq, start=0, step=5, window=20, k=3):
    """ Takes a sequence and calculates k-mer probabilities in it.
    Then, splits the sequence in a series of sliding windows
    and calculates the k-mer based Shannon entropy using those probabilities.
    It is different from calculating the probabilities for each individual window.
    

    :param seq: A sequence represented as a string.
    :type seq: str
    :param start: Position to start the sliding windows (0-indexed).
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
    Kmers_pi=dict()
    Kmers_values=dict()

    # going through the sequence in forward direction
    for w_start in range(start, len(seq) - k + 1):
        kmer = seq[w_start:w_start+k]

        if kmer in Kmers:
            Kmers[kmer] = Kmers[kmer] + 1
        else:
            Kmers[kmer] = 1
    
    # making reverse complement
    sequence=Seq(seq)
    sequence_revcomp=sequence.reverse_complement()
    seq_revcomp=str(sequence_revcomp)

    # going through the sequence in reverse direction
    for w_start in range(0, len(seq) - k - start + 1):
        kmer = seq_revcomp[w_start:w_start+k]

        if kmer in Kmers:
            Kmers[kmer] = Kmers[kmer] + 1
        else:
            Kmers[kmer] = 1


    # calculating pi for each kmer
    total_instances=np.sum(np.fromiter(Kmers.values(), dtype=float))
    for key in Kmers:
        kmer_pi=float(Kmers[key])/total_instances
        Kmers_pi[key]=kmer_pi
        Kmers_values[key]=kmer_pi*np.log2(kmer_pi)

    # calculating entropy for each window, using pi values calculated for the entire DNA sequence, both strands
    for w_start in range(start, len(seq) - window + 1, step):
        w_end = min(w_start + window, len(seq))

        curr_kmers=[]

        for start in range(w_start, w_end - k + 1):
            kmer = seq[start:start+k]
            if kmer not in curr_kmers:
                curr_kmers.append(kmer)

        window_entropy=-np.sum(Kmers_values[x] for x in curr_kmers)

        Hk_seq.append([w_start, w_end, window_entropy])

    return Hk_seq

