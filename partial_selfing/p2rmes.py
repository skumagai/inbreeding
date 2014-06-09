# -*- mode: python; coding: utf-8; -*-

# p2rmes.py - brief description

# Copyright (C) 2014 Seiji Kumagai

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

from __future__ import print_function

import argparse as ap, sys
from os.path import split

def genotype(l):
    g = l.split('\t')[1:]
    return ' '.join(str(int(g[i] != g[i+1])) for i in range(0, len(g), 2))


def write(ifname, ofname):
    with open(ifname, 'r') as ifh, open(ofname, 'w') as ofh:
        # force Windows newline '\r\n' instead of Unix' '\n'
        ofh.write("1\r\n")

        # put population headers first
        n = split(ifh.name)[1]
        s, l = [s.strip() for s in [next(ifh), next(ifh)]]
        next(ifh)

        ofh.write('{}\r\n{}\r\n{}\r\n'.format(n, s, l))

        for l in ifh:
            ofh.write(genotype(l.strip()) + '\r\n')


if __name__ == '__main__':

    p = ap.ArgumentParser(description='convert phase-format genotype data to RMES input format')
    p.add_argument('INFILE',
                   type=str)
    args = p.parse_args()
    ofname = split(args.INFILE)[1][:-3] + 'rmes.txt'
    write(args.INFILE, ofname)
