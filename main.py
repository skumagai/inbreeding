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
            float(t['value']) / (4 * N)
            for t in self._p['mutation']['theta']
            for rep in range(t['times'])
            ]

        try:
            self._pp = cobj['post process']
        except:
            pass

        self.gens = self._p['N'] * self._g['gens']
        self.burnin = self._p['N'] * self._g['burnin']
        self.output_per = self._p['N'] * self._g['output per']

        self.mode = self._p['mutation']['model']
        mating = self._p['mating']
        self.model = mating['model']
        try:
            self.a = float(mating['a'])
            self.tau = float(mating['tau'])
            self.sigma = float(mating['sigma'])
        except:
            pass

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
    config = Config(json.load(args.CONFIG), args.STR)
    import partial_selfing.post_process as pp
    pp.run(config)


if __name__ == '__main__':

    cp = argparse.ArgumentParser(add_help=False)
    cp.add_argument("CONFIG",
                    type=argparse.FileType('r'),
                    help='configuration of simulation')
    cp.add_argument('STR',
                    type=str,
                    nargs='*',
                    help='string substituted in output file name')

    parser = argparse.ArgumentParser()

    subparsers = parser.add_subparsers()
    s = subparsers.add_parser('simulate',
                              help='run simulations',
                              parents=[cp])
    s.set_defaults(func=simulate)

    pp = subparsers.add_parser('post-process',
                               help='post process simulated data',
                               parents=[cp])

    pp.set_defaults(func=post_process)

    args = parser.parse_args()
    args.func(args)
