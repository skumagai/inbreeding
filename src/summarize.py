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

def parseArgs()
    parser = argparse.ArgumentParser(
        description='compute f, g, P, W')

    parser.add_argument('DIR',
                        nargs='*',
                        help='name of directories storing simuPOP.Population object')
    return parser.parse_args()


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


def get_info(dir):
    with open(os.path.join(dir, 'conf.json'), 'r') as rf:
        info = json.load(rf)
    return info


def chunks(l, n):
    '''Divide a list into equal-length sublists'''
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def compute_f(pop, nloci):
    nsites = pop.totNumLoci() / nloci
    fs = [0] * nloci
    pop_size = pop.popSize()
    for ind in pop.individuals():
        geno0 = chunks(ind.genotype(ploidy=0), nsites)
        geno1 = chunks(ind.genotype(ploidy=1), nsites)
        fs = check_identity(fs, geno0, geno1)
    return [float(val) / pop_size for val in fs]


def comput_g(pop, nloci):
    nsites = pop.totNumLoci() / nloci
    gs = [0] * nloci
    pop_size = pop.popSize()
    for pair in itertools.combinations(pop.individuals(), 2):
        geno00 = chunks(pair[0].genotype(ploidy=0), nsites)
        geno01 = chunks(pair[0].genotype(ploidy=1), nsites)
        geno10 = chunks(pair[1].genotype(ploidy=0), nsites)
        geno11 = chunks(pair[1].genotype(ploidy=1), nsites)
        gs = check_identity(gs, geno00, geno10)
        gs = check_identity(gs, geno00, geno11)
        gs = check_identity(gs, geno01, geno10)
        gs = check_identity(gs, geno01, geno11)
    denom = 4. * pop_size * (pop_size - 1) / 2
    return [float(val) / denom for val in gs]


def compute_P(pop, nloci):
    nsites = pop.totNumLoci() / nloci
    Ps = {key: 0 for key in itertools.product(xrange(2), repeat = nloci)}
    zeros = [0] * nloci
    pop_size = pop.popSize()
    for ind in pop.individuals():
        geno0 = chunks(ind.genotype(ploidy=0), nsites)
        geno1 = chunks(ind.genotype(ploidy=1), nsites)
        Ps[tuple(check_identity(zeros, geno0, geno1))] += 1
    return {key: float(val) / pop_size for key, val in Ps.iteritems()}


def compute_W():
    nsites = pop.totNumLoci() / nloci
    Ws = {key: 0 for key in itertools.product(xrange(2), repeat = nloci)}
    zeros = [0] * nloci
    pop_size = pop.popSize()
    for pair in itertools.combinations(pop.individuals(), 2):
        geno00 = chunks(pair[0].genotype(ploidy=0), nsites)
        geno01 = chunks(pair[0].genotype(ploidy=1), nsites)
        geno10 = chunks(pair[1].genotype(ploidy=0), nsites)
        geno11 = chunks(pair[1].genotype(ploidy=1), nsites)
        Ws[tuple(check_identity(zeros, geno00, geno10))] += 1
        Ws[tuple(check_identity(zeros, geno00, geno11))] += 1
        Ws[tuple(check_identity(zeros, geno01, geno10))] += 1
        Ws[tuple(check_identity(zeros, geno01, geno11))] += 1
    denom = 4. * pop_size * (pop_size - 1) / 2
    return {key: float(val) / denom for key, val in Ws.iteritems()}


def check_identity(counts, geno0, geno1):
    for i, geno_pair in enumerate(zip(geno0, geno1)):
        if geno_pair[0] == geno_pair[1]:
            counts[i] += 1
    return counts


def summarise(d, mode):
    info = get_info(d)
    num_loci = info[u'number of loci']
    for f in [str(s) for s in info[u'files']]:
        pop = sim.loadPopulation(f)
        f = compute_f(pop, num_loci)
        g = compute_g(pop, num_loci)
        P = compute_P(pop, num_loci)
        W = compute_W(pop, num_loci)


def main():
    args = parseArgs()
    mode = import_right_module(args)
    for d in args.DIR:
        f, g, P, W = summarise(d, mode)
    sys.exit(0)

if __name__ == '__main__':
    main()
