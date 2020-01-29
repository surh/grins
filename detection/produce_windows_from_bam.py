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


import pysam
import argparse


class windows:
    """A list of window elements with methods to ensure that new overlapping
    windows are merged into existing ones."""

    def __init__(self):
        self.windows = []

    def find_overlaps(self, w):
        """Finds if any existing window overlaps with window w, and
        returns list of indexes of overlapping windows"""

        overlaps = []
        for i in range(len(self.windows)):
            check1 = w.start < self.windows[i].end
            check2 = self.windows[i].start < w.end
            if check1 and check2:
                overlaps.append(i)

        return overlaps

    def update_window(self, i, start, end):
        """Update existing window.
        Probably should be a method of class window"""
        self.windows[i].start = start
        self.windows[i].end = end

    def remove_windows(self, index):
        """Remove a list of windows by index"""

        # Important to reverse to ensure the correct windows are removed
        for i in sorted(index, reverse=True):
            del self.windows[i]

    def add_window(self, w):
        """Adds a new window. If no overlaps are found, just add the window.
        If overlaps. Define new window by merging all overlaps, and then
        replace old windows with the new merged one."""

        overlaps = self.find_overlaps(w=w)
        if len(overlaps) == 0:
            self.windows.append(w)
        elif len(overlaps) > 0:
            start = w.start
            end = w.end
            for i in overlaps:
                start = min(start, self.windows[i].start)
                end = max(end, self.windows[i].end)

            new_w = window(start, end)
            self.remove_windows(overlaps)
            self.add_window(w=new_w)
        else:
            raise ValueError("Incorrect number of overlaps")

    def n_windows(self):
        """Returns number of windows in the list"""

        return len(self.windows)


class window:
    """A simple window class"""

    def __init__(self, start, end):
        if start < 0 or end < 0:
            print("Start: ", start, "End: ", end)
            raise ValueError("start and end must be non negative.")
        if type(start) is not int or type(end) is not int:
            raise ValueError("start and end must be integers.")
        if start > end:
            raise ValueError("start cannot be greater than end.")
        self.start = start
        self.end = end


def process_arguments():
    # Read arguments
    parser_format = argparse.ArgumentDefaultsHelpFormatter
    parser = argparse.ArgumentParser(formatter_class=parser_format)
    required = parser.add_argument_group("Required arguments")

    # Define description
    parser.description = ("Find windows of repeated elements from bam file")

    # Define required arguments
    required.add_argument("--input", help=("A bam file as produced from "
                                           "converting bowtie2 output. "
                                           "Query names must be START_END."),
                          required=True, type=str)

    # Define other arguments
    parser.add_argument("--output", help=("Name for output GFF3 file"),
                        type=str,
                        default="grins.gff3")
    parser.add_argument("--w_size", help=("Window size (needed?)"),
                        type=int,
                        default=150)

    # Read arguments
    print("Reading arguments")
    args = parser.parse_args()

    # Processing goes here if needed

    return args


def find_bam_windows(file, min_size=150):
    """Reads bam file and finds alignment windows"""
    # Read samfile and produce list of alignment windows
    samfile = pysam.AlignmentFile(file, "rb")
    reads = samfile.fetch(until_eof=True)
    multi_windows = []
    for r in reads:
        q_ref, q_start, q_end = r.query_name.split('|')
        q_start = int(q_start)
        q_end = int(q_end)

        r_ref, r_start, r_end = r.reference_name.split('|')
        r_start = int(r_start)
        r_end = int(r_end)

        aln_len = min(q_end - q_start, r_end - r_start)

        # Keep minimum aln length and depending on record
        # discard self maps.
        if aln_len >= min_size and r.pos > 0:
            if q_ref != r_ref:
                my_read = [q_ref, q_start, q_end,
                           r_ref, r_start, r_end,
                           r.mapping_quality]
                multi_windows.append(my_read)
            else:
                if r.pos != q_start:
                    my_read = [q_ref, q_start, q_end,
                               r_ref, r_start, r_end,
                               r.mapping_quality]
                    multi_windows.append(my_read)



    return multi_windows


def merge_bam_windows(bam_windows, w_size=150):
    """Needs improvement. Takes output from find_bam_windows and produces
    windows object with merged overlapping windows.
    NOTE: No support for multiple references (i.e. contig/chromosome)
    """

    res = windows()
    res = dict()
    for w in bam_windows:
        # Add query window
        if w[0] not in res:
            res[w[0]] = windows()
        new_w = window(start=w[1], end=w[2])
        res[w[0]].add_window(w=new_w)

        # add target
        if w[3] not in res:
            res[w[3]] = windows()
        new_w = window(start=w[4], end=w[5])
        res[w[3]].add_window(w=new_w)

    return res


def write_gff3(windows, output):
    with open(output, 'w') as oh:
        oh.write("##gff-version 3\n")
        i = 1
        for k in windows:
            for w in windows[k].windows:
                id = ''.join(['ID=pGRINS_', str(i)])
                # Create row, switch to 1-indexed with closed interval
                row = [k, 'grinspred', 'pGRINS',
                       str(w.start+1), str(w.end),
                       '.', '+', '.', id]
                row = "\t".join(row)
                oh.write(row + "\n")
                i = i + 1
    oh.close()

    return output


if __name__ == "__main__":
    args = process_arguments()
    bam_windows = find_bam_windows(args.input, min_size=args.w_size)
    # print(bam_windows)
    res_windows = merge_bam_windows(bam_windows=bam_windows,
                                    w_size=args.w_size)
    print("Found:", res_windows.n_windows(), "windows")
    write_gff3(windows=res_windows, output=args.output)
