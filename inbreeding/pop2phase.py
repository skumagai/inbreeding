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
import os.path
import operator

from inbreeding import utility

# From https://github.com/onyxfish/csvkit/commit/5f2c5d9d5c8596a4c7440e53dba0bd8d92feb3b9
# Ensure SIGPIPE doesn't throw an exception
# Prevents [Errno 32] Broken pipe errors, e.g. when piping to 'head'
try:
    import signal
    signal.signal(signal.SIGPIPE, signal.SIG_DFL)
except ImportError:
    pass


### Functions specific to 'csv' subcommand.
# In order to choose between 'long' and 'binary' modules, for
# respectively the infinite-alleles and infinite-sites models, command
# line arguments have to be processed before the import.
def parseArgs():
    parser = argparse.ArgumentParser(
        description="convert simuPOP's Population object to genotype data format used in phase.")

    subparsers = parser.add_subparsers(help='sub-command help')
    csv_p = subparsers.add_parser('csv', help='generate CSV file of genotypes')
    csv_p.add_argument('DIR',
                       type=str,
                       help='directory')
    csv_p.add_argument('-o', '--output',
                       metavar='FILE',
                       type=argparse.FileType('w'),
                       default=sys.stdout,
                       help='save converted data into FILE')
    csv_p.set_defaults(func=do_csv)


    phase_p = subparsers.add_parser('phase', help='generate phase file of genotypes')
    phase_p.add_argument('SAMPLE',
                         type=int,
                         help='number of individuals to sample')
    phase_p.add_argument('LOCI',
                         type=str,
                         help='number of loci to sample (comma separated list)')
    phase_p.add_argument('CSVS',
                         type=str,
                         help='name of CSV files storing genotypes (comma separaeted list)')
    phase_p.add_argument('-r', '--random',
                         action='store_true',
                         help='randomize order of samples')
    phase_p.add_argument('-o', '--output',
                         metavar='FILE',
                         type=argparse.FileType('w'),
                         default=sys.stdout,
                         help='save converted data into FILE')
    phase_p.set_defaults(func=do_phase)

    return parser.parse_args()


def get_mode(args):
    info = utility.get_info(args.DIR)

    mode = 'infinite-sites'
    if 'mode' in info and info['mode'] == 'infinite-alleles':
        mode = 'infinite-alleles'
    return mode


def encode(locus, codes):
    try:
        length = codes['length']
    except KeyError:
        codes['length'] = length = len(locus)

    if length != len(locus):
        raise ValueError('wrong number of sites in a locus: exp {}, obs{}'.
                         format(codes['length'], length))
    locus = str(locus)

    try:
        code = codes[locus]
    except KeyError:
        codes[locus] = code = codes['idx']
        codes['idx'] += 1
    return code


def convert_genotype(genotype, loci_dict):
    return [encode(locus, loci_dict[idx])
            for idx, locus in enumerate(genotype)]


def to_csv(args, mode):
    d = args.DIR
    info = utility.get_info(d)

    # this is needed as simuPOP is loaded in a different file.  It's
    # name is not bound in the namescope of this file.
    sim = sys.modules['simuPOP']

    pop = sim.loadPopulation(os.path.join(d, str(info['files'][0])))
    pop_size = pop.popSize()
    nloci = info['number of loci']
    lenlocus = pop.totNumLoci() / nloci

    loci_dict = {i:{'idx': 0} for i in range(pop_size)}
    genotypes = []

    for ind in pop.individuals():
        geno0 = convert_genotype(utility.chunks(ind.genotype(ploidy = 0), lenlocus), loci_dict)
        geno1 = convert_genotype(utility.chunks(ind.genotype(ploidy = 1), lenlocus), loci_dict)

        genotypes.append(itertools.chain.from_iterable(list(zip(geno0, geno1))))
    return genotypes


def write_csv(genotypes, output):
    pop_size = len(genotypes)
    writer = csv.writer(output)
    writer.writerow([pop_size])

    for genotype in genotypes:
        writer.writerow(list(genotype))


### Functions specific to 'phase' subcommand.
def get_csvs(args):
    files = [os.path.expanduser(f) for f in args.CSVS.split(',')]
    for f in files:
        if os.path.isfile(f) is False:
            print('[ERROR] {} is not a file.'.format(f), file=sys.stderr)
            sys.exit(1)
    return files


def to_phase(args):
    csvs = get_csvs(args)
    sloci = [int(i) for i in args.LOCI.split(',')]
    if len(sloci) != len(csvs):
        print('[ERROR] inconsistent sampling strategy.', file=sys.stderr)
        sys.exit(1)

    inds = [[] for i in range(args.SAMPLE)]
    for c, l in zip(csvs, sloci):
        ind_data = sample_loci(c, l, args.SAMPLE)
        for i, d in enumerate(ind_data):
            inds[i].extend(d)
    return inds


def sample_loci(c, nloci, nsample):
    try:
        with open(c, 'r') as rf:
            reader = csv.reader(rf)
            pop_size = int(reader.next()[0])
            try:
                idx = random.sample(range(pop_size), nsample)
            except (ValueError, TypeError) as e:
                print('[ERROR] {}'.format(e), file=sys.stderr)
                sys.exit(1)
            inds = [l for i, l in enumerate(reader) if i in idx]
    except IOError as e:
        print('[ERROR] {}'.format(e), file=sys.stderr)
        sys.exit(1)

    totalloci = len(inds[0]) / 2
    order = [(i * 2, i * 2 + 1) for i in random.sample(range(totalloci), nloci)]
    extractor = operator.itemgetter(*sorted(itertools.chain.from_iterable(order)))
    loci = [extractor(ind) for ind in inds]
    return loci


def write_phase(nsample, nloci, data, output):
    nloci = sum([int(i) for i in nloci.split(',')])

    print(nsample, file=output)
    print(nloci, file=output)
    print('M' * nloci, file=output)
    for i, d in enumerate(data):
        print("sample_{}\t".format(i), file=output, end='')
        print("\t".join([str(di) for di in d]), file=output)


def randomize(data):
    nloci = len(data[0])
    order = range(nloci)
    random.shuffle(order)
    randomizer = operator.itemgetter(*order)
    genotypes = [randomizer(ind) for ind in data]
    return genotypes


### (Almost) top-level functions
def do_csv(args):
    '''convert genotypic information in simuPOP.Population object into CSV file.'''
    mode = utility.import_right_module(args)
    data = to_csv(args, mode)
    write_csv(data, args.output)
    sys.exit(0)


def do_phase(args):
    '''convert genotypic information generated by pop2phase.py csv [...] into phase file.'''
    data = to_phase(args)
    if args.random is True:
        data = randomize(data)
    write_phase(args.SAMPLE, args.LOCI, data, args.output)
    sys.exit(0)


def main():
    args = parseArgs()
    args.func(args)


if __name__ == '__main__':
    main()
