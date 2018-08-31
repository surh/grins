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

# from Bio.Seq import Seq
from Bio.SeqFeature import FeatureLocation
from functools import singledispatch


def calculate_skew(seq, which='GC'):
    """Calculate GC or AT skew from Seq or string object"""

    # Count nucleotides
    count_a = seq.upper().count('A')
    count_c = seq.upper().count('C')
    count_g = seq.upper().count('G')
    count_t = seq.upper().count('T')

    # Make sure there are no ambiguous bases
    base_count = count_a + count_c + count_g + count_t
    if base_count != len(seq):
        raise ValueError("Sequence base count differe from sequence length!")

    # Calculate skew
    if which == 'GC':
        skew = (count_g - count_c) / (count_g + count_c)
    elif which == 'AT':
        skew = (count_a - count_t) / (count_a + count_t)
    elif which == 'CG':
        skew = (count_c - count_g) / (count_g + count_c)
    elif which == 'TA':
        skew = (count_t - count_a) / (count_a + count_t)
    else:
        raise ValueError("Invalid skew type ({})".format(which))

    return skew


@singledispatch
def calculate_window_skew(location, seq, which='GC'):
    """Calculate skew in a window defined by location.
    Generic function default"""

    # Get sequence and make sure it falls in range
    my_window = location.extract(seq)
    if len(location) != len(my_window):
        raise IndexError("The location given falls out of range")

    try:
        skew = calculate_skew(my_window, which=which)
    except ValueError:
        raise
    except ZeroDivisionError:
        raise

    return skew


@calculate_window_skew.register(list)
@calculate_window_skew.register(tuple)
def _(location, seq, which='GC'):
    """Calculate skew in a window defined by two coordinates.
    List implementation, requires a list of two integers."""

    # Check list/tuple size
    if len(location) != 2:
        raise ValueError("The location list should have exactly two elements")

    # Create feature location object and call default generic method
    loc = FeatureLocation(location[0], location[1])
    try:
        skew = calculate_window_skew(loc, seq, which)
    except ValueError:
        raise
    except ZeroDivisionError:
        raise

    return skew


def calculate_sliding_skew(seq, w_size=100, w_step=1, which='GC'):
    """Calculate skew with across a sliding window"""

    # Iterate over windows
    Tab = []
    for i in range(0, len(seq) - w_size + 1, w_step):
        # Define FeatureLocation object and call generic method
        loc = FeatureLocation(i, i + w_size)
        skew = calculate_window_skew(loc, seq=seq, which=which)
        Tab.append([loc, skew])

    return Tab
