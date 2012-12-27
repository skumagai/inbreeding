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

from __future__ import print_function

import argparse
import csv
import itertools
import json
import sys
import random

# in order to choose between 'long' and 'binary' modules, for
# respectively the infinite-alleles and infinite-sites models, command
# line arguments have to be processed before the import.
def parseArgs():
    parser = argparse.ArgumentParser(
        description="convert simuPOP's Population object to genotype data format used in phase.")

    subparsers = parser.add_subparsers()

    # command-line arguments for storing genotype data of an entire
    # population into CSV file.
    prepare_parser = subparsers.add_parser('prepare',
                                           help='save genotype into CSV')
    prepare_parser.add_argument('DIR',
                                nargs='*',
                                type=str,
                                help='name of directories storing simuPOP.Population object')
    prepare_parser.set_defaults(func=prepare)

    conv_parser = subparsers.add_parser('generate',
                                        help='generate phase files')

    conv_parser.add_argument('-r', '--random',
                             action='store_true',
                             help='randomize order of samples')
    conv_parser.add_argument('-s', '--sample-size',
                             type=int,
                             nargs='?',
                             help='sample size from each CSV file')
    conv_parser.add_argument('-o', '--output',
                             metavar='FILE',
                             type=argparse.FileType('w'),
                             default=sys.stdout,
                             help='save converted data into FILE')
    conv_parser.add_argument('-i', '--input',
                             metavar='FILE',
                             nargs='*',
                             type=argparse.FileType('r'),
                             default=sys.stdin,
                             help='list of input files')
    conv_parser.set_defaults(func=generate)

    return parser.parse_args()


def get_info(dir):
    with open(os.path.join(dir, 'conf.json'), 'r') as rf:
        info = json.load(rf)
    return info

def extract_data_rep(dir, info):
    num_loci = info[u'number of loci']
    data = []
    for pop in [sim.loadPopulation(os.path.join(dir, str(filename)))
                for filename in info[u'files']]:
        num_sites = pop.totalNumLoci() / num_loci
        if pop.totalNumLoci() % num_loci != 0:
            print('[ERROR] total number of sites is not multiple of ' +
                  'sites per locus {}'.format(num_loci),
                  file=sys.stderr)
            sys.exit(1)
        locus_dict = [{'idx': 1} for i in xrange(num_loci)]
        genotypes = [[] for i in xrange(num_loci)]

        for ind in pop.individuals():
            geno0 = ind.genotype(ploidy = 0)
            geno1 = ind.genotype(ploidy = 1)
            for locus in xrange(num_loci):
                start = locus * num_sites
                stop = (locus + 1) * num_sites
                l0_code = encode(geno0[start:stop], locus_dict[locus])
                l1_code = encode(geno1[start:stop], locus_dict[locus])
                genotypes[locus].append((l0_code, l1_code))

    return data


def encode(locus, codes):
    locus = str(locus)
    try:
        code = codes[locus]
    except:
        codes[locus] = codes['idx']
        code = codes['idx']
        codes['idx'] += 1
    return codes




# entry point to subcommand 'prepare'
def prepare(args):
    import an appropriate
    try:
        import simuOpt
        info = get_info(args.DIR[0])
        if u'mode' in info and info[u'mode'] == 'infinite-allele':
            simuOpt.setOptions(alleleType = 'long')
        else:
            simuOpt.setOptions(alleleType = 'binary')
        import simuPOP as sim
    except:
        print('[ERROR] prepare command needs simuPOP', file=sys.stderr)
        sys.exit(1)

    for directory in args.DIR:
        data = get_data(directory)
        try:
            values[triplet].append(data)
        except:
            values[triplet] = [data]

        with open(directory + '.csv', 'w') as wf:
            writer = csv.writer(wf)
            writer.writerow(['population size', len(data[0])])
            writer.writerow(['mutation rate', triplet[0]])
            writer.writerow(['recombination rate', triplet[1]])
            wirter.writerow(['selfing rate', triplet[2]])
            for d in zip(*data):
                writer.writerow(itertools.chain(*d))


# entry point to subcommand 'generate'
def generate(args):



def main():
    args = parseArgs()
    args.func(args)
    sys.exit(0)

    # if args:
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

if __name__ == '__main__':
    main()
