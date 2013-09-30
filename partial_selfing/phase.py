# -*- mode: python; coding: utf-8; -*-

# phase.py - output PHASE format haplotype data from simulation results.

# Copyright (C) 2013 Seiji Kumagai

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


import argparse, csv, random

def process(infile, outfile, nsam):
    with open(infile, 'r') as cfile:
        reader = csv.reader(cfile)

        # header
        reader.next()

        idx = 1

        row = reader.next()
        npop = int(row[1])
        nloci = int(row[4])

        samples = random.sample(xrange(2 * npop), nsam)

        with open(outfile, 'w') as ocfile:
            ocfile.write("{}\n{}\n".format(nsam, nloci))
            ocfile.write("{}\n".format('M'*nloci))

            if idx in samples:
                ocfile.write("sample_{}\t{}\n".format(idx, '\t'.join(row[-nloci:])))

            for row in reader:
                idx += 1
                if idx in samples:
                    ocfile.write("sample_{}\t{}\n".format(idx, '\t'.join(row[-nloci:])))



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('infile')
    parser.add_argument('outfile')
    parser.add_argument('nsam')

    args = parser.parse_args()

    process(args.infile, args.outfile, int(args.nsam))
