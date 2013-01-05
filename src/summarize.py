# -*- mode: python; coding: utf-8; -*-

# summarize.py - brief description

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
import json
import csv
import os.path
import sys
import itertools

from utility import import_right_module, chunks, get_info

def parseArgs():
    parser = argparse.ArgumentParser(
        description='compute f, g, P, W')

    parser.add_argument('DIR',
                        nargs='*',
                        help='name of directories storing simuPOP.Population object')
    parser.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='output file')

    return parser.parse_args()


def basic_info(inds):
    return len(inds), len(inds[0][0])


def add_lists(lst0, lst1):
    if len(lst0) != len(lst1):
        raise IndexError('Length mismatch')
    for i, v in enumerate(lst1):
        lst0[i] += v


def compute_f(inds):
    pop_size, nloci = basic_info(inds)
    fs = [0] * nloci
    for ind in inds:
        add_lists(fs, check_identity(nloci, ind[0], ind[1]))
    return [float(val) / pop_size for val in fs]


def compute_g(inds):
    pop_size, nloci = basic_info(inds)
    gs = [0] * nloci
    for pair in itertools.combinations(inds, 2):
        geno00, geno01 = pair[0]
        geno10, geno11 = pair[1]
        add_lists(gs, check_identity(nloci, geno00, geno10))
        add_lists(gs, check_identity(nloci, geno00, geno11))
        add_lists(gs, check_identity(nloci, geno01, geno10))
        add_lists(gs, check_identity(nloci, geno01, geno11))
    denom = 4. * pop_size * (pop_size - 1) / 2
    return [float(val) / denom for val in gs]


def compute_P(inds):
    pop_size, nloci = basic_info(inds)
    Ps = {key: 0 for key in itertools.product(xrange(2), repeat = nloci)}
    for ind in inds:
        Ps[tuple(check_identity(nloci, ind[0], ind[1]))] += 1
    return {key: float(val) / pop_size for key, val in Ps.iteritems()}


def compute_W(inds):
    pop_size, nloci = basic_info(inds)
    Ws = {key: 0 for key in generate_PW_keys(nloci)}
    for pair in itertools.combinations(inds, 2):
        geno00, geno01 = pair[0]
        geno10, geno11 = pair[1]
        Ws[tuple(check_identity(nloci, geno00, geno10))] += 1
        Ws[tuple(check_identity(nloci, geno00, geno11))] += 1
        Ws[tuple(check_identity(nloci, geno01, geno10))] += 1
        Ws[tuple(check_identity(nloci, geno01, geno11))] += 1
    denom = 4. * pop_size * (pop_size - 1) / 2
    return {key: float(val) / denom for key, val in Ws.iteritems()}


def check_identity(nloci, geno0, geno1):
    counts = [0] * nloci
    for i, geno_pair in enumerate(zip(geno0, geno1)):
        if geno_pair[0] == geno_pair[1]:
            counts[i] = 1
    return counts


def generate_PW_keys(nloci):
    return itertools.product(xrange(2), repeat = nloci)


def summarise(d, mode):
    info = get_info(d)
    num_loci = info[u'number of loci']
    keys = list(generate_PW_keys(num_loci))

    # this is needed as simuPOP is loaded in a different file.  It's
    # name is not bound in the namescope of this file.
    sim = sys.modules['simuPOP']

    mut= info[u'mutation rate'][u'scaled']
    rec = info[u'recombination rate'][u'scaled']
    selfing = info[u'selfing rate']

    results = []
    header = ['mutation rate', 'recombination rate', 'selfing rate', 'replicate id']
    header += ['f_' + str(i) for i in xrange(num_loci)]
    header += ['g_' + str(i) for i in xrange(num_loci)]
    header += ['P_' +  ''.join(str(i) for i in xrange(num_loci))
               + '(' + ''.join(str(k) for k in key) + ')' for key in keys]
    header += ['W_' +  ''.join(str(i) for i in xrange(num_loci))
               + '(' + ''.join(str(k) for k in key) + ')' for key in keys]

    for i, f in enumerate([str(s) for s in info[u'files']]):
        pop = sim.loadPopulation(os.path.join(d, f))
        pop_size = pop.popSize()
        nsites = pop.totNumLoci() / num_loci

        inds = [[tuple(chunks(ind.genotype(ploidy=0), nsites)),
                 tuple(chunks(ind.genotype(ploidy=1), nsites))] for ind in pop.individuals()]

        f = compute_f(inds)
        g = compute_g(inds)
        P = compute_P(inds)
        W = compute_W(inds)
        results.append([mut, rec, selfing, i] +
                       f + g + [P[key] for key in keys] + [W[key] for key in keys])

    return header, results


def main():
    args = parseArgs()
    mode = import_right_module(args)

    writer = csv.writer(args.output)
    write_header = True
    for d in args.DIR:
        header, results = summarise(d, mode)
        if write_header is True:
            nelems = len(results)
            writer.writerow(header)
            write_header = False
        elif len(results) != nelems:
            print('[WARN] write header in the middle of CSV file.', file=sys.stderr)
            writer.writerow(header)

        for result in results:
            writer.writerow(result)

    sys.exit(0)

if __name__ == '__main__':
    main()
