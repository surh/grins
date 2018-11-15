#!/usr/bin/env python
# Copyright (C) 2018 Sur Herrera Paredes

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


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Read tab-delimited file with column and row "
                          "names, and rename them")

    # Define required arguments
    required.add_argument("--infile", help=("Input file"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--outfile", help=("Output file"),
                        type=str,
                        default="outfile.txt")
    parser.add_argument("--mapfile", help=("ID map file"),
                        type=str,
                        default="mapfile.txt")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def create_col_names_dict(col_names):
    """Create dictionary from list"""

    NEWIDS = dict()
    newids = []
    i = 1
    for name in col_names:
        new_name = 'G' + "%9d" % i
        if name in NEWIDS:
            raise ValueError("Repeated ID")
        NEWIDS[name] = new_name
        newids.append(new_name)
        i = i + 1

    print(newids)
    return NEWIDS, newids


if __name__ == "__main__":
    args = process_arguments()

    with open(args.infile, 'r') as ih, open(args.outfile, 'w') as oh:
        # Read, process and write header
        print("Processing header")
        header = ih.readline()
        col_names = header.split("\t")
        NEWIDS, new_header = create_col_names_dict(col_names)
        oh.write("\t".join([new_header]) + "\n")

        # Process following lines
        print("Processing lines")
        for line in ih:
            fields = line.split("\t")
            if fields[0] in NEWIDS:
                fields[0] = NEWIDS[fields[0]]
            else:
                raise ValueError("Row ID not found ({})".format(fields[0]))

            new_line = "\t".join([fields]) + "\n"
            oh.write(new_line)
    ih.close()
    oh.close()

    # Write map file
    print("Writing map file")
    with open(args.mapfile, 'w') as oh:
        for k in NEWIDS:
            line = "\t".join([k, NEWIDS[k]]) + "\n"
            oh.write(line)
    oh.close()
