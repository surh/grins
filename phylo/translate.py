#!/usr/bin/env python

import argparse
from Bio import SeqIO
from Bio.Alphabet import generic_dna
import warnings


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Takes a file with sequence records and translates "
                          "each one of them.")

    # Define required arguments
    required.add_argument("--infile", help=("Name of the input file with DNA "
                                            "sequence records"),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--check_cds", help=("If passed, the program will "
                                             "check that each record is a "
                                             "full CDS, meaning it starts "
                                             "with a valid start codon "
                                             "(translated to methionine), "
                                             "that the sequence is a "
                                             "multiple of three, and that "
                                             "there is a valid stop codon "
                                             "at the end"),
                        action="store_true",
                        default=False)
    parser.add_argument("--dont_stop", help=("By default, the translation "
                                             "will stop at any stop codon. "
                                             "Passing this argument will "
                                             "continue pass any stop codon."),
                        action="store_true",
                        default=False)
    parser.add_argument("--format", help=("Format of input and output files"),
                        type=str, default="fasta")
    parser.add_argument("--outfile", help=("Name of the output file with "
                                           "translated sequence records"),
                        type=str,
                        default="translated.faa")
    parser.add_argument("--remove_stops", help=("If passed, any trailing stop "
                                                "codons will be removed from "
                                                "the translated sequence"),
                        action="store_true",
                        default=False)
    parser.add_argument("--rename", help=("Rename sequences so they have "
                                          "short IDs"),
                        action="store_true",
                        default=False)
    parser.add_argument("--rename_table", help=("Mapping file with the new "
                                                "sequence names."),
                        default="sequence_names_map.txt",
                        type=str)
    parser.add_argument("--table", help=("Genetic code table, must be a "
                                         "name or NCBI identifier"),
                        type=str,
                        default="Standard")

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed
    args.to_stop = not args.dont_stop

    if(args.format != 'fasta'):
        warnings.warn(("WARNING: The script has only been tested with fasta "
                       "sequences. While use of other formats is supported "
                       "by Bio.SeqIO, it might lead to unexpected results."))

    return args


if __name__ == "__main__":
    args = process_arguments()

    NAMES = dict()
    with open(args.infile, 'r') as ih, open(args.outfile, 'w') as oh:

        i = 0
        for record in SeqIO.parse(ih, args.format,
                                  alphabet=generic_dna):
            prot = record.seq.translate(table=args.table,
                                        to_stop=args.to_stop,
                                        cds=args.check_cds)

            if args.remove_stops:
                prot = prot.rstrip("*")

            # Change record
            record.seq = prot

            # Count after translation
            i = i + 1

            # Rename if needed
            if args.rename:
                new_name = 'S' + "%07d" % i
                NAMES[new_name] = [record.id, record.description]
                record.id = new_name
                record.description = new_name

            # Write
            SeqIO.write(record, oh, args.format)
    ih.close()
    oh.close()
    print("Processed {} sequence records.".format(str(i)))
    if args.rename:
        with open(args.rename_table, 'w') as ih:
            print("Writing {}".format(args.rename_table))
            i = 0
            for id in NAMES:
                line = "\t".join([id, NAMES[id][0], NAMES[id][1]]) + "\n"
                ih.write(line)
                i = i + 1
        ih.close()
        print("Wrote {} sequence names.".format(str(i)))
