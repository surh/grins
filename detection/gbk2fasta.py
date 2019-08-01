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
from Bio import SeqRecord
import argparse
import os


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Convert gbk to fasta.")

    # Define required arguments
    required.add_argument("--input", help=("Input filename in genbank"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Output filename. If empty, the "
                                          "output will be "
                                          "<basename>.fasta"),
                        type=str,
                        default="")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()
    i = 0
    for record in SeqIO.parse(args.input, args.format):
        if i > 0:
            raise ValueError("There should be only one record in your input. "
                             "Only the first record was processed.")
        print(record.id)

        # Create output name if needed
        if args.output == '':
            name = os.path.basename(os.path.splitext('/home/sur/test.txt')[0])
            args.output = name + '.fasta'

        with open(args.output, 'w') as oh:
            SeqIO.write(record, oh, "fasta")
        oh.close()
        i = i + 1
