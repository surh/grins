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

import argparse
import os
import shutil


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Script that runs throuhg a directory structure "
                          "containing genomes, and creates a parallel "
                          "structure in a specified folder, discarding "
                          "unwanted genomes")

    # Define required arguments
    required.add_argument("indir", help=("Directory with genome direcetories"),
                          type=str)

    # Define other arguments
    parser.add_argument("--keep", help=("File with genomes to keep, one per "
                                        "line"),
                        type=str)
    parser.add_argument("--discard", help=("File with genomes to discard, "
                                           "one per line"),
                        type=str)
    parser.add_argument("--outdir", help=("Directory where to place cleaned "
                                          "genomes."),
                        type=str,
                        default="output/")
    parser.add_argument("--records", help=("Filename of action records"),
                        type=str, default='records.txt')
    parser.add_argument("--copy", help=("Copy wanted genomes instead of link"),
                        default=False, action="store_true")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    if args.keep is None and args.discard is None:
        raise ValueError("At least one of keep and discard must be provided")

    return args


def list_level_dirs(path, level=1):
    """Recursive function to return all directories at a specified level"""

    if level == 1:
        # Base case, when current level is wanted
        dirs = os.listdir(path)
    elif level > 1:
        # If we need to go deeper
        next_level_dirs = os.listdir(path)
        dirs = []
        for d in next_level_dirs:
            tempd = list_level_dirs(os.path.join(path, d), level - 1)
            dirs.extend([os.path.join(d, t) for t in tempd])
    else:
        # Error
        raise ValueError("ERROR: level must be a positive integer")

    return dirs


def read_list(path):
    """Read file into list, one value per line"""

    with open(path, 'r') as ih:
        res = []
        for line in ih:
            res.append(line.rstrip("\n"))
    ih.close()

    return(res)


def keep_dir(src_dir, dest_dir, copy=False):
    """Make symlink or copy of dir"""
    # Check if the dest_dir already exists

    if not os.path.isdir(dest_dir):
        os.mkdir(dest_dir)

    dirname = os.path.basename(src_dir)
    if copy:
        shutil.copytree(src_dir, os.path.join(dest_dir, dirname))
    else:
        os.symlink(os.path.realpath(src_dir),
                   os.path.join(dest_dir, dirname))

    return


if __name__ == "__main__":
    args = process_arguments()

    keep = []
    discard = []
    print("Reading list of genomes to keep/discard.")
    if args.keep is not None:
        keep = read_list(args.keep)
    if args.discard is not None:
        discard = read_list(args.discard)

    print("Reading genome dirs.")
    genome_dirs = list_level_dirs(args.indir, 2)

    print("Preparing output directory.")
    if os.path.isdir(args.outdir):
        raise ValueError("The output directory already exists.")
    else:
        os.mkdir(args.outdir)

    print("Processing genome dirs")
    Record = dict()
    for gdir in genome_dirs:
        genome = os.path.basename(gdir)
        print("\tEvaluating genome {} at {}".format(genome, gdir))

        group_dir = os.path.dirname(gdir)
        src_dir = os.path.join(args.indir, gdir)
        dest_dir = os.path.join(args.outdir, group_dir)

        keep_flag = genome in keep
        discard_flag = genome in discard

        # Consider all the cases
        if len(keep) > 0 and len(discard) > 0:
            # Both keep and discard specifications
            if keep_flag and discard_flag:
                # genome marked for opposite actions
                raise ValueError("Genome {} is marked for both "
                                 "discard and keep")
            elif keep_flag:
                print("\t\tKeeping...")
                keep_dir(src_dir, dest_dir, args.copy)
                Record[genome] = [gdir, group_dir, 'Kept']
            elif discard_flag:
                print("\t\tDiscarding...")
                Record[genome] = [gdir, group_dir, 'Discarded']
            else:
                # Equivalent to discard but with different record
                raise ValueError("Genome not included in keep or discard")
                Record[genome] = [gdir, group_dir, 'Skipped']
        elif len(keep) > 0:
            # Only keep specification
            if keep_flag:
                print("\t\tKeeping...")
                keep_dir(src_dir, dest_dir, args.copy)
                Record[genome] = [gdir, group_dir, 'Kept']
            else:
                print("\t\tDiscarding...")
                Record[genome] = [gdir, group_dir, 'Discarded']
        elif len(discard) > 0:
            # Only discard specification
            if discard_flag:
                print("\t\tDiscarding...")
                Record[genome] = [gdir, group_dir, 'Discarded']
            else:
                print("\t\tKeeping...")
                keep_dir(src_dir, dest_dir, args.copy)
                Record[genome] = [gdir, group_dir, 'Kept']
        else:
            # No specification
            raise ValueError("At least one of keep or discard must be passed")

    print("Writing records")
    with open(args.records, 'w') as oh:
        for k in Record:
            fields = Record[k]
            fields = [k] + fields
            fields = "\t".join(fields) + "\n"
            oh.write(fields)
    oh.close()
