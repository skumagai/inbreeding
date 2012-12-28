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

# in order to choose between 'long' and 'binary' modules, for
# respectively the infinite-alleles and infinite-sites models, command
# line arguments have to be processed before the import.
def parseArgs():
    parser = argparse.ArgumentParser(
        description="convert simuPOP's Population object to genotype data format used in phase.")

    parser.add_argument('SAMPLE',
                        type=int,
                        help='sample size')
    parser.add_argument('LOCI',
                        type=str,
                        help='number of samples (command separated list)')
    parser.add_argument('DIR',
                        nargs='*',
                        type=str,
                        help='name of directories storing simuPOP.Population object')
    parser.add_argument('-r', '--random',
                             action='store_true',
                             help='randomize order of samples')
    parser.add_argument('-o', '--output',
                        metavar='FILE',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='save converted data into FILE')
    return parser.parse_args()


def get_info(dir):
    with open(os.path.join(dir, 'conf.json'), 'r') as rf:
        info = json.load(rf)
    return info


def encode(locus, codes):
    locus = str(locus)
    try:
        code = codes[locus]
    except:
        codes[locus] = codes['idx']
        code = codes['idx']
        codes['idx'] += 1
    return code


def import_right_module(args):
    mode = 'infinite-sites'
    try:
        global sim
        import simuOpt
        info = get_info(args.DIR[0])
        if u'mode' in info and info[u'mode'] == u'infinite-alleles':
            mode = 'infinite-alleles'
            simuOpt.setOptions(alleleType = 'long')
        else:
            simuOpt.setOptions(alleleType = 'binary')
        import simuPOP as sim
    except:
        print('[ERROR] needs simuPOP', file=sys.stderr)
        sys.exit(1)
    return mode

def check_directories_exist(dirs, mode):
    ds = []
    for d in dirs:
        if os.path.isdir(d) is True:
            with open(os.path.join(d, 'conf.json')) as rf:
                if json.load(rf)[u'mode'] == mode:
                    ds.append(d)
                else:
                    print('[ERROR] cannot mix mutation models.', file=sys.stderr)
                    sys.exit(1)
        else:
            print('[ERROR] {} is not a directory.'.format(d), file=sys.stderr)
            sys.exit(1)
    return ds


def convert(args, mode):
    dirs = check_directories_exist(args.DIR, mode)
    loci = [int(i) for i in args.LOCI.split(',')]
    if len(loci) != len(dirs):
        print('[ERROR] inconsistent sampling strategy.', file=sys.stderr)
        sys.exit(1)

    inds = [[] for i in xrange(args.SAMPLE)]
    for d, l in zip(dirs, loci):
        ind_data = sample_loci(d, l, args.SAMPLE)
        for i, d in enumerate(ind_data):
            inds[i].extend(d)
    print(inds)
    return inds

def sample_loci(d, nloci, nsample):
    info = get_info(d)
    info = get_info(d)
    pop_size = info[u'population size']
    num_loci = info[u'number of loci']
    pop = sim.loadPopulation(os.path.join(d, str(info[u'files'][0])))
    allele_len = pop.totNumLoci() / num_loci

    ind_idx = random.sample(xrange(pop_size), nsample)
    inds = [pop.individual(idx) for idx in ind_idx]
    loci_idx = random.sample(xrange(num_loci), nloci)
    loci_dict = {i:{'idx': 0} for i in loci_idx}

    genotypes = []

    for ind in inds:
        geno0 = chunks(ind.genotype(ploidy = 0), allele_len)
        geno0 = [encode(locus, loci_dict[idx])
                 for idx, locus in enumerate(geno0)
                 if idx in loci_idx]
        geno1 = chunks(ind.genotype(ploidy = 1), allele_len)
        geno1 = [encode(locus, loci_dict[idx])
                 for idx, locus in  enumerate(geno1)
                 if idx in loci_idx]

        genotypes.append((geno0, geno1))

    genotypes = [list(itertools.chain(*[list(i) for i in zip(*ind)]))
                 for ind in genotypes]
    return genotypes


def chunks(l, n):
    '''Divide a list into equal-length sublists'''
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def write_data(nsample, nloci, data, output):
    nloci = sum([int(i) for i in nloci.split(',')])

    print(nsample, file=output)
    print(nloci, file=output)
    print('M' * nloci, file=output)
    for i, d in enumerate(data):
        print("sample_{}\t".format(i), file=output, end='')
        print("\t".join([str(di) for di in d]), file=output)


def main():
    args = parseArgs()
    mode = import_right_module(args)
    data = convert(args, mode)
    write_data(args.SAMPLE, args.LOCI, data, args.output)
    sys.exit(0)


if __name__ == '__main__':
    main()
