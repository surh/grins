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

# Load libraries
from Bio import SeqIO
from Bio import SeqRecord
import argparse
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Split sequence into fasta file windows.")

    # Define required arguments
    required.add_argument("--input", help=("Input filename"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Output filename. If empty, the "
                                          "output will be "
                                          "<basename>_windows.fasta"),
                        type=str,
                        default="")
    parser.add_argument("--w_size", help=("Window size"),
                        type=int,
                        default=150)
    parser.add_argument("--s_size", help=("Step size for the sliding window"),
                        type=int,
                        default=30)
    parser.add_argument("--format", help=("Format of input file"),
                        type=str,
                        default='genbank',
                        choices=['fasta', 'genbank'])

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def create_window_records(record, w_size=150, s_size=30):
    """Create window records and return as list of records. Record
    IDs match start and end positions"""
    windows = []
    r_len = len(record)
    for pos in range(0, r_len, s_size):
        end = min(pos+w_size, r_len)
        window = record[pos:end]
        w_id = '_'.join([str(pos), str(end)])
        windows.append(SeqRecord.SeqRecord(id=w_id, seq=window.seq,
                                           description=''))

    return windows


if __name__ == "__main__":
    args = process_arguments()
    i = 0
    for record in SeqIO.parse(args.input, args.format):
        if i > 0:
            raise ValueError("There should be only one record in your input. "
                             "Only the first record was processed.")
        print(record.id)
        windows = create_window_records(record, args.w_size, args.s_size)

        # Create output name if needed
        if args.output == '':
            name = os.path.basename(os.path.splitext(args.input)[0])
            args.output = name + '_windows.fasta'

        with open(args.output, 'w') as oh:
            SeqIO.write(windows, oh, "fasta")
        oh.close()
        i = i + 1
