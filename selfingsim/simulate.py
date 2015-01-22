"""
selfingsim.simulate
###################

Run a forward-in-time population genetics simulation for partially selfing organisms.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# standard imports
import argparse
import json
import sys

def run():
    """
    Runs simulations as a stand-alone script.
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    setup_command_line(subparsers)
    args = parser.parse_args()
    args.func(args)

def setup_command_line(subparsers):
    """
    Sets up command line interface.
    """
    # Forward simulation
    parser = subparsers.add_parser(
        'simulate',
        help='run forward simulations with selfing')
    parser.add_argument(
        'config',
        type=argparse.FileType('r'),
        help='the settings of simulation as a JSON file')
    parser.add_argument(
        'result',
        type=str,
        help='the name of file storing simulation results')
    parser.set_defaults(func=simulate)

def simulate(args):
    """
    Runs simulation.
    """
    config = Config(json.load(args.CONFIG), args.STR)

    if config.mode == 'infinite sites':
        exec_infinite_sites(config)
    elif config.mode == 'infinite alleles':
        exec_infinite_alleles(config)
    else:
        print("Unknown mutational model specified", file=sys.stderr)
        sys.exit(1)

class Config(object):
    """
    Stores input parameters for simulations include mutation rates, population
    size etc.
    """

    def __init__(self, cobj, subst):
        N = int(cobj['population']['N'])
        self._g = cobj['general']
        self._p = cobj['population']
        self.outfile = self._g['outfile'].format(*subst)
        self.m = [
            float(t['value']) / (4 * N)
            for t in self._p['mutation']['theta']
            for _ in range(t['times'])]

        try:
            self._pp = cobj['post process']
        except KeyError:
            pass

        self.gens = self._p['N'] * self._g['gens']
        self.burnin = self._p['N'] * self._g['burnin']
        self.output_per = int(self._p['N'] * self._g['output per'])

        self.mode = self._p['mutation']['model']
        mating = self._p['mating']
        self.model = mating['model']

        try:
            self.s = float(mating['s'])
        except NameError:
            self.a = float(mating['a'])
            self.tau = float(mating['tau'])
            self.sigma = float(mating['sigma'])
            at = self.a * self.tau

            if self.model == 'androdioecy':
                self.s = at / (at + (1 - self.a) * self.sigma)
            else:
                Nh = N * float(self._p['sex ratio'])
                Nf = N - Nh
                self.s = at * Nh / (at * Nh + Nh * (1 - self.a) \
                        + Nf * self.sigma)
                self.h = Nh * (1 - self.a) / (Nh * (1 - self.a) \
                        + Nf * self.sigma)

        try:
            self.allele_length = self._p['allele_length']
        except KeyError:
            self.allele_length = 1


    def __getattr__(self, name):
        name1 = name.replace('_', ' ')
        try:
            return self._g[name1]
        except KeyError:
            try:
                return self._p[name1]
            except KeyError:
                try:
                    return self._pp[name1]
                except KeyError:
                    raise KeyError('{}'.format(name))


# Selectively import simuPOP with an appropriate alleleType.  Because
# python caches imported modules, second or later call of import
# does not simuPOP re-imported.  Therefore, those
# submodules will automatically use the right version of simuPOP.
def exec_infinite_sites(config):
    """
    Launches simulations with the infinite-sites model.
    """
    import selfingsim.infinite_sites as model
    model.run(config)


# See the comment in front of exec_two_loci
def exec_infinite_alleles(config):
    """
    Launches simulations with the infinite-alleles model.
    """
    import selfingsim.infinite_alleles as model
    model.run(config)

if __name__ == '__main__':
    run()
