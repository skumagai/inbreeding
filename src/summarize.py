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

    subparsers = parser.add_subparsers(help='sub-command help')

    stat_p = subparsers.add_parser('stats', help='compute f, g, P, and W')


    stat_p.add_argument('DIR',
                        nargs='*',
                        help='name of directories storing simuPOP.Population object')
    stat_p.add_argument('-o', '--output',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='output file')

    stat_p.set_defaults(func=stats)

    plot_p = subparsers.add_parser('plot', help='plot f, g, P, and W')

    plot_p.add_argument('CSV',
                        type=str,
                        help='input CSV file')
    plot_p.set_defaults(func=plot)


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
    # Eventually multiple {f,g}^{*} columns are merged when computing
    # mean and sem.
    header += ['f^{' + str(i) + '}' for i in xrange(num_loci)]
    header += ['g^{' + str(i) + '}'  for i in xrange(num_loci)]
    header += ['P_{' + str(num_loci) + '}('
               + ''.join(str(k) for k in key) + ')' for key in keys]
    header += ['W_{' + str(num_loci) + '}('
               + ''.join(str(k) for k in key) + ')' for key in keys]

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


def stats(args):
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


def filter_column_labels(data, cond):
    return [key for key in data.columns if cond(key)]


def mathtext(fn):
    def _fn(*args):
        return r'$' + fn(*args) + '$'
    return _fn


@mathtext
def get_ylabel(label):
    first_letter = label[0][0]
    if first_letter == 'f' or  first_letter == 'g':
        return first_letter
    else:
        idx = label.find('(')
        return label[:idx]

@mathtext
def get_label(label):
    return label


def plot_category(ax, index, group, keys, func):
    xs = filter_column_labels(group, lambda x: (x[0] in keys and x[1] == 'mean'))
    xerrs = filter_column_labels(group, lambda x: (x[0] in keys and x[1] == 'sem'))
    colors = ['b', 'g', 'r', 'c', 'm', 'y', 'k', 'w']
    for x, xerr, col in zip(xs, xerrs, colors):
        label = get_label(x[0])
        ax.errorbar(index, group[x], yerr=group[xerr] / 2, fmt='-', color=col, label=label)
        ax.plot(index, [func(s) for s in index], 'd', color=col)
    ax.set_ylabel(get_ylabel(keys[0]))
    ax.legend()


def plot(args):
    global pd, np, scipy, plt
    import pandas as pd
    import numpy as np
    import scipy.stats
    import matplotlib.pyplot as plt

    summary_funcs = [np.mean, scipy.stats.sem]

    raw_data = pd.read_csv(args.CSV, header=0, index_col=list(xrange(4)))
    raw_data = raw_data.sort_index()
    fs = filter_column_labels(raw_data, lambda x: x[0] == 'f')
    gs = filter_column_labels(raw_data, lambda x: x[0] == 'g')
    Ps = filter_column_labels(raw_data, lambda x: x[0] == 'P')
    Ws = filter_column_labels(raw_data, lambda x: x[0] == 'W')

    for (mut, rec), group in \
      raw_data.groupby(level=['mutation rate', 'recombination rate']):

      group = group.reset_index(level=['mutation rate', 'recombination rate'],
                                drop=True)

      merged_f = pd.concat(group[key] for key in fs).groupby(level='selfing rate').agg(summary_funcs)
      merged_f.rename(columns={'mean': ('f', 'mean'), 'sem': ('f', 'sem')}, inplace=True)

      merged_g = pd.concat(group[g] for g in gs).groupby(level='selfing rate').agg(summary_funcs)
      merged_g.rename(columns={'mean': ('g', 'mean'), 'sem': ('g', 'sem')}, inplace=True)

      for keys in (fs, gs):
        for key in keys:
          del group[key]

      group = group.groupby(level='selfing rate').agg(summary_funcs)
      group = pd.concat([group, merged_f, merged_g], axis=1)

      index = group.index
      func = None

      # prepare 2x2 panel for plotting f, g, P, and W.
      fig, ((ax_f, ax_g), (ax_P, ax_W)) = plt.subplots(2, 2, sharex='all')

      # urgly hack to share x-axis label
      invisible_ax = fig.add_subplot(111)
      invisible_ax.set_frame_on(False)
      invisible_ax.set_xticks([])
      invisible_ax.set_yticks([])
      invisible_ax.set_xlabel('selfing rate', labelpad=20)

      # plot f
      f_func = lambda s: 1 - (1 - s) * mut / (1 + (1 - s / 2) * mut)
      plot_category(ax_f, index, group, ['f'], f_func)

      # # plot g
      g_func = lambda s: 1 - (1 - s / 2) * mut / (1 + (1 - s / 2) * mut)
      plot_category(ax_g, index, group, ['g'], g_func)

      # # plot P
      # P_func = lambda
      # plot_category(ax_P, index, group, Ps, None)

      # # plot W
      # plot_category(ax_W, index, group, Ws, None)

      fig.suptitle(r'$\theta = {}, r = {}$'.format(mut, rec))

    plt.show()


def main():
    args = parseArgs()
    args.func(args)

if __name__ == '__main__':
    main()
