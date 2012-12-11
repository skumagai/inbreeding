# -*- mode: python; coding: utf-8; -*-

# pop2x.py - convert simuPOP's Population object to something else

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
import sys
import random

# in order to choose between 'long' and 'binary' modules, for
# respectively the infinite-alleles and infinite-sites models, command
# line arguments have to be processed before the import.
def parseArgs():
    parser = argparse.ArgumentParser(description="convert simuPOP's Population object to a usable form to something else.")
    parser.add_argument('SIZE',
                        type=int,
                        help='sample size')
    parser.add_argument('LOCI',
                        type=int,
                        help='number of sampled loci')
    parser.add_argument('INPUT',
                        type=str,
                        help='input file containing Population object (.pop)')
    parser.add_argument('OUTPUT',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='output file')
    parser.add_argument('--num-segre-sites',
                        metavar='NSITES',
                        type=int,
                        default=256,
                        help='maximum number of segregating sites per locus (default: 256)')
    parser.add_argument('--infinite-alleles',
                        action='store_true',
                        help='use the infinite-alleles model instead of the infinite-sites model')
    return parser.parse_args()

args = parseArgs()

import simuOpt
if args.infinite_alleles == True:
    mode = 'long'
else:
    mode = 'binary'

simuOpt.setOptions(alleleType = mode)

import simuPOP as sim

if __name__ == '__main__':
    pop = sim.loadPopulation(args.INPUT)
    output = args.OUTPUT

    if args.infinite_alleles == True:
        allele_len = 1
    else:
        allele_len = args.num_segre_sites

    tot_num_loci = pop.totNumLoci()
    allele = args.num_segre_sites

    if allele_len * args.LOCI > tot_num_loci or tot_num_loci % allele_len != 0:
        sys.exit('inconsistent chromosomal length')

    output.write(str(args.SIZE) + '\n')
    output.write(str(args.LOCI) + '\n')
    output.write(str('M' * args.LOCI + '\n'))

    sampled_loci = random.sample(xrange(tot_num_loci / allele_len), args.LOCI)

    loci_dict = {i: {'index': 1} for i in sampled_loci}

    for i, ind in enumerate([pop.individual(i)
                             for i in random.sample(xrange(pop.popSize()), args.SIZE)]):
        output.write('sample {}'.format(i))
        for locus in sampled_loci:
            for ploidy in xrange(2):
                l = str(ind.genotype(ploidy = ploidy)[
                    (locus * allele_len):((locus + 1) * allele_len)])
                try:
                    output.write('\t' + str(loci_dict[locus][l]))
                except KeyError:
                    loci_dict[locus][l] = loci_dict[locus]['index']
                    loci_dict[locus]['index'] += 1
                    output.write('\t' + str(loci_dict[locus][l]))
        output.write('\n')
