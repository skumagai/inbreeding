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


def write(a):
    fhs = a.INFILE
    # force Windows newline '\r\n' instead of Unix' '\n'
    print(len(fhs), end='\r\n')

    # put population headers first
    for fh in fhs:
        n = split(fh.name)[1]
        s, l = [s.strip() for s in [next(fh), next(fh)]]
        next(fh)

        print('{}\r\n{}\r\n{}\r\n'.format(n, s, l), end='')

    for fh in fhs:
        for l in fh:
            print(genotype(l.strip()), end='\r\n')


if __name__ == '__main__':

    p = ap.ArgumentParser(description='convert phase-format genotype data to RMES input format')
    p.add_argument('INFILE',
                   nargs='*',
                   default=[sys.stdin],
                   type=ap.FileType('r'))
    write(p.parse_args())
