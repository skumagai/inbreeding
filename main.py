# -*- mode: python; coding: utf-8; -*-

# main.py - An entry point to partial selfing simulations.  The
# simulations are conducted with population genetics forward simulator
# simuPOP, and four scenarios are explored.

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

import argparse, json, sys


class Config(object):

    def __init__(self, cobj, subst):
        N = cobj['population']['N']
        self._g = cobj['general']
        self._p = cobj['population']
        self.outfile = self._g['outfile'].format(*subst)
        self.m = [
            float(t['value']) / N
            for t in self._p['mutation']['theta']
            for rep in range(t['times'])
            ]

        try:
            self._pp = cobj['post process']
        except:
            pass

        self.gens = self._g['N'] * self._g['gens']
        self.burnin = self._g['N'] * self._g['burnin']

        self.mode = self._p['mutation']['model']
        self.model = self._p['mating']['model']

        m = self._p['mating']
        if self.model == 'pure hermaphrodite':
            self.s = m['a'] * m['tau']
            self.s = self.s / (self.s + 1 - m['a'])
        elif self.model == 'androdioecy':
            self.s = self.s / (self.s + (1 - m['a']) * m['sigma'])
        else:
            Nh = N * self._p['sex ratio']
            Nf = N - Nh
            sg = m['tau'] * Nh * m['a']
            sg = sg / (sg + Nh * (1 - m['a']) + Nf * m['sigma'])
            h = Nh * (1 - m['a'])
            h = h / (h + Nf * m['sigma'])
            self.s = (sg, h)

        try:
            self.allele_length = self._p['allele_length']
        except:
            self.allele_length = 1


    def __getattr__(self, name):
        name1 = name.replace('_', ' ')
        try:
            return self._g[name1]
        except:
            try:
                return self._p[name1]
            except:
                try:
                    return self._pp[name1]
                except:
                    raise KeyError('{}'.format(name))



# Selectively import simuPOP with an appropriate alleleType.  Because
# python caches imported modules, subsequent import of simuPOP in
# model specific submodules is not re-imported.  Therefore, those
# submodules will automatically use the right version of simuPOP.
def exec_infinite_sites(config):
    import partial_selfing.infinite_sites as model
    model.run(config)


# See the comment in front of exec_two_loci
def exec_infinite_alleles(config):
    import partial_selfing.infinite_alleles as model
    model.run(config)


def summarize(config):
    import partial_selfing.summary as model
    model.run(config)

def simulate(args):
    config = Config(json.load(args.CONFIG), args.STR)

    if config.mode == 'infinite sites':
        exec_infinite_sites(config)
    else:
        exec_infinite_alleles(config)

def post_process(args):
    pass

if __name__ == '__main__':

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()
    s = subparsers.add_parser('simulate',
                              help='run simulations')
    s.set_defaults(func=simulate)
    s.add_argument("CONFIG",
                   type=argparse.FileType('r'),
                   help='configuration of simulation')
    s.add_argument('STR',
                   type=str,
                   nargs='*',
                   help='string substituted in output file name')
    pp = subparsers.add_parser('post-process',
                               help='post process simulation data')
    # pp.set_defaults(func=post_process)


    args = parser.parse_args()
    args.func(args)
