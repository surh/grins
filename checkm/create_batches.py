#!/usr/bin/env python
# Copyright (C) 2018-2019 Sur Herrera Paredes

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


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("It splits a collection of genomes into batches of "
                          "the specified size. "
                          "It creates a directory structure at the specified "
                          "outdir with one subdirectory per batch and "
                          "symlinks within each subdirectory to absolute "
                          "paths of the "
                          "genome files in the corresponding batch. "
                          "It writes a mapping file that groups genomes "
                          "into batches.\n"
                          "It assumes that the directory passed has a set "
                          "of subdirectories corresponding to some arbitrary "
                          "group. Then each group subdirectory has one "
                          "subdirectory per genome where the only .fna "
                          "file in the genome subdirectory corresponds to "
                          "the contigs file.")

    # Define required arguments
    required.add_argument("--indir", help=("Directory that containes genomes. "
                                           "See above for the expected "
                                           "directory structure."),
                          required=True, type=str)
    required.add_argument("--outdir", help=("Directory where to create "
                                            "symbolic link structure"),
                          type=str,
                          required=True)

    # Optional arguments
    parser.add_argument("--outfile", help=("Name of the mapping file to be "
                                           "created."),
                        type=str, default='batch_map.txt')
    parser.add_argument("--batch_size", help=("Number of genomes per batch"),
                        required=True, type=int, default=100)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    return args


def find_genome_files(indir):
    """Read indir and find all genome files. Check assumptions
    about tree structure"""

    Genomes = []
    specname_dirs = os.listdir(indir)
    for spec in specname_dirs:
        if os.path.isdir(os.path.join(indir, spec)):
            print("\tReading genome dirs from {}".format(spec))
            genome_dirs = os.listdir(os.path.join(indir, spec))
            for genome in genome_dirs:
                genome_path = os.path.join(indir,
                                           spec,
                                           genome,
                                           genome + '.fna')
                if os.path.isfile(genome_path):
                    Genomes.append(genome_path)
                else:
                    raise FileNotFoundError("Genome file ({}) not "
                                            "found".format(genome_path))
        else:
            print("{} is not a directory. Skipping.")

    return Genomes


def make_batches(Genomes, outdir, batch_size):
    """Make batches"""

    # Create output directory
    os.mkdir(outdir)

    # Iterate over genomes
    curr_batch = 0
    curr_batch_size = 0
    curr_batch_dir = ''
    batch_name = ''
    Map = []
    for genome in Genomes:
        if (curr_batch_size % batch_size) == 0:
            # Initialize batch dir
            batch_name = ''.join(['batch_', str(curr_batch)])
            curr_batch_dir = ''.join([outdir, '/', batch_name, '/'])
            os.mkdir(curr_batch_dir)
            curr_batch = curr_batch + 1
            curr_batch_size = 0
            print("\tProcessing batch {}".format(batch_name))

        # Create symlink
        # print(genome, curr_batch_dir)
        src_abspath = os.path.abspath(genome)
        genome_name = os.path.basename(genome)
        target_name = ''.join([curr_batch_dir, '/', genome_name])
        target_name = os.path.abspath(target_name)
        os.symlink(src=src_abspath, dst=target_name)
        Map.append([batch_name, src_abspath, target_name])
        curr_batch_size = curr_batch_size + 1
        # print("======", curr_batch_size, batch_size)

    return Map


if __name__ == "__main__":
    args = process_arguments()

    # Get list of genome files
    print("Getting list of genome files")
    try:
        Genomes = find_genome_files(args.indir)
        # print(Genomes)
    except:
        raise

    print("Making batches")
    Map = make_batches(Genomes, args.outdir, args.batch_size)

    print("Writting mapping file")
    # print(Map)
    with open(args.outfile, 'w') as oh:
        for l in Map:
            line = '\t'.join(l)
            line = ''.join([line, "\n"])
            oh.write(line)
    oh.close()
