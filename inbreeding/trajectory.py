# -*- mode: python; coding: utf-8; -*-

# trajectory.py - aim to help deciding how long burn-in period should
# be by plotting time series

# Copyright (C) 2012 Seiji Kumagai

# Permission is hereby granted, free of charge, to any person obtaining a
# copy of this software and associated documentation files (the "Software"),
# to deal in the Software without restriction, including without limitation
# the rights to use, copy, modify, merge, publish, distribute, sublicense,
# and/or sell copies of the Software, and to permit persons to whom the
# Software is furnished to do so, subject to the following conditions:

# The above copyright notice and this permission notice (including the next
# paragraph) shall be included in all copies or substantial portions of the
# Software.

# THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
# IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
# FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT.  IN NO EVENT SHALL
# THE AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
# LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
# FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
# DEALINGS IN THE SOFTWARE.

import argparse
import matplotlib.pyplot as plt
import pandas
import sys


def parse_args():
    parser = argparse.ArgumentParser(
        description='plot heterozygosity or number of segregating sites over generations')
    parser.add_argument('INPUT',
                        type=argparse.FileType('r'),
                        help='time series data in CSV format')
    parser.add_argument('-c', '--columns',
                        metavar='RANGE',
                        action='append',
                        type=str,
                        help='use RANGE of columns to plot')
    return parser.parse_args()

def get_range(columns, mincol, maxcol):
    ranges = []
    if len(columns) == 0:
        return [(mincol, maxcol)]
    for col in columns:
        start, stop = col.split('-')
        start = int(start)
        stop = int(stop)
        if mincol <= start < maxcol and mincol < stop <= maxcol:
            ranges.append((start, stop))
        else:
            sys.exit('index out of range')
    return ranges

if __name__ == '__main__':
    args = parse_args()

    data = pandas.read_csv(args.INPUT)

    nmin = 2
    nmax = data.shape[1]

    for start, stop in get_range(args.columns, nmin, nmax):
        data.ix[:,start:(stop-1)].plot()
        lgd = plt.legend()
        lgd.set_visible(False)
    plt.show()
