#!/usr/bin/env python
# Copyright (C) 2020 Sur Herrera Paredes

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
import json


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes output from antiSMASH 5 and creates a "
                          "GFF3 file of all the regions.")

    # Define required arguments
    required.add_argument("--input",
                          help=("JSON file produced by antiSMASH 5."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output",
                        help=("Filename for output GFF3 file. If empty then "
                              "the input file name will be used replacing the "
                              ".json extension for .gff3"),
                        type=str,
                        default="")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


if __name__ == "__main__":
    args = process_arguments()

    # Read json
    with open(args.input, 'r') as ih:
        asmash = json.load(ih)
