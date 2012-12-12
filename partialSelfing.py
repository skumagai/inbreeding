# -*- mode: python; coding: utf-8; -*-

# inbreeding.py - simulation of two-linked loci under the infinite-site model

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
from collections import deque
import sys
import random
import os.path

# Command line arguments have to be processed before importing
# simuOpt, because the right module to import depends on what mutation
# model is used.  Furthermore, a mutation model is specified via a
# command line flag.  Without any additional flag, the infinite-sites
# model is assumed.  With --infinite-alleles flag, the
# infinite-alleles model is used.
#
# In the former, 'binary' version of simuPOP is required, and 'long'
# version of simuPOP is required in case of the latter.
def parseArgs():
    parser = argparse.ArgumentParser(description="run partial selfing simulations")
    parser.add_argument('POP',
                        type=int,
                        help='population size')
    parser.add_argument('NGEN',
                        type=int,
                        help='number of generations excluding burn-in period')
    parser.add_argument('NREP',
                        type=int,
                        help='number of replicates')
    parser.add_argument('OUTPUT',
                        type=str,
                        help='save final populations')
    parser.add_argument('-r', '--recombination-rate',
                        metavar='RHO',
                        type=float,
                        help='scaled recombination rate (default: free recombinaton)')
    parser.add_argument('-s', '--selfing-rate',
                        metavar='S',
                        type=float,
                        default=0,
                        help='selfing rate (default: 0)')
    parser.add_argument('-m', '--mutation-rate',
                        metavar=('POS', 'THETA'),
                        type=str,
                        nargs=2,
                        action='append',
                        help='positions and scaled rates of mutation. (default:["all", "0"])')
    parser.add_argument('-b', '--burn-in',
                        metavar='B',
                        type=int,
                        default=0,
                        help='burn-in (default: 0)')
    parser.add_argument('-n', '--num-loci',
                        metavar='NLOCI',
                        type=int,
                        default=2,
                        help='number of loci (default: 2)')
    parser.add_argument('--num-segre-sites',
                        metavar='NSITES',
                        type=int,
                        default=256,
                        help='maximum number of segregating sites per locus (default: 256)')
    parser.add_argument('--seed',
                        type=int,
                        default=0,
                        help='random number seed (default: use posix time)')
    parser.add_argument('--trajectory',
                        metavar='TRAJ_FILE',
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=None,
                        const=sys.stdout,
                        help='record heterozygosity and number of segregating sites each generation for later inspection to determine an appropriate durtion of burn-in.  (default: STDOUT)')
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

simuOpt.setOptions(alleleType = mode)# ,
                   # optimized = False,
                   # quiet = True)

import simuPOP as sim

# Unify the handling of result reporting.
class Writer(sim.PyOperator):
    '''Interface of writing output'''

    def __init__(self, output, num_loci, allele_len, burnin, *args, **kwargs):
        self.num_loci = num_loci
        self.allele_len = allele_len
        self.output = output
        self.burnin = burnin
        self.has_header_printed = False
        super(Writer, self).__init__(func = self.write, *args, **kwargs)

    def write(self, pop):
        '''stub for a function to write data'''
        pass

    def write_header(self):
        '''stub for a function to write header'''
        pass

    def _write_common(self, pop):
        dvars = pop.dvars()
        self.output.write('{},{},'.format(dvars.rep, dvars.gen - self.burnin))
        self.output.write(','.join(['{}'.format(val)
                                    for val in computeHeterozygosity(pop,
                                                                     self.num_loci,
                                                                     self.allele_len)]))

    def _write_common_header(self):
        if self.has_header_printed == False:
            self.output.write('"rep","gen",')
            self.output.write(','.join(['"1-f (locus {})"'.format(i)
                                        for i in xrange(self.num_loci)]))
            self.has_header_printed = True



class InfSiteWriter(Writer):
    '''Write heterozigosities and number of segregating sites at each generation.

    This class is intended to be used in explorative runs and only
    under the infinite-sites model.'''

    def __init__(self, output, num_loci, allele_len, burnin, *args, **kwargs):
        super(InfSiteWriter, self).__init__(output,
                                            num_loci,
                                            allele_len,
                                            burnin,
                                            *args,
                                            **kwargs)

    def write(self, pop):
        self._write_common(pop)
        self.output.write(','.join(['{}' for i in xrange(self.num_loci)]).format(
            *computeNumberOfSegregatingSites(pop, num_loci, allele_len)) + '\n')
        return True

    def write_header(self):
        self._write_common_header()
        self.output.write(',' + ','.join(['"# segre (locus {})"'.format(i)
                                          for i in xrange(self.num_loci)]) + '\n')


class InfAlleleWriter(Writer):
    '''Write heterozigosities at each generation.

    This class is intended to be used in explorative runs and only
    under the infinite-alleles model.'''

    def __init__(self, output, num_loci, burnin, *args, **kwargs):
        super(InfAlleleWriter, self).__init__(output,
                                              num_loci,
                                              1,
                                              burnin,
                                              *args,
                                              **kwargs)

    def write(self, pop):
        self._write_common(pop)
        self.output.write('\n')
        return True

    def write_header(self):
        self._write_common_header()
        self.output.write('\n')


# Implement infinite-sites mutation model.  As new mutations arise,
# number of polymorphic sites increases.  In the current
# implementation a fix number of slots are pre-allocated to store
# polymorphic sites.  Once all available space is exhausted, sites
# that are monomorphic (due to random drift) are recycled.  When the
# recycle is failed to open up any new space, simulations are
# terminated.

#
# Dynamically increasing the number of slots is planned but not
# implemented.

class Mutator(sim.PyOperator):
    '''base class for custom Mutator classes'''
    def __init__(self, mu, num_loci, allele_len, rep, burnin, *args, **kwargs):
        self.num_loci = num_loci
        self._build_mutation_rates(mu)
        self.allele_len = allele_len
        self.burnin = burnin
        super(Mutator, self).__init__(func = self.mutate, *args, **kwargs)

    def mutate(self, pop):
        '''stub for mutation method.  To be implemented in children.'''
        raise NotImplementedError


    def _build_mutation_rates(self, mut):
        '''build per-locus mutation rate usable in Mutator classes.

        Rates can be specified for all loci (all), arange of loci (e.g. 0-10),
        or list of loci (e.g. 0,2,4).  If more than one entries are present for
        a single locus, the last one is effective.'''

        num_loci = self.num_loci
        mut_rates = [0.0] * num_loci

        for m in mut:
            loc, val = m
            val = float(val)
            parts = loc.split(',')
            for part in parts:
                if part == 'all':
                    mut_rates = [val] * num_loci
                elif part.find('-') >= 0:
                    boundaries = part.split('-')
                    if len(boundaries) != 2:
                        sys.exit('loc spec range must has the form "start-stop"')
                    start, stop = [int(d) for d  in boundaries]
                    # stop is not inclusive following python's convension.
                    if 0 <= start < num_loci and 0 <= stop <= num_loci:
                        for pos in xrange(start, stop):
                            mut_rates[pos] = val
                    else:
                        sys.exit('mutation loc spec: range out of loci')
                else:
                    pos = int(part)
                    if 0 <= pos < num_loci:
                        mut_rates[int(part)] = val
                    else:
                        sys.exit('mutation loc spec: position out of loci')
        self.mu = mut_rates


class InfSiteMutator(Mutator):

    def __init__(self, mu, num_loci, allele_len, rep, burnin, *args, **kwargs):
        self.available = {r: [deque(xrange(allele_len))
                              for i in range(num_loci)]
                          for r in range(rep)}
        super(InfSiteMutator, self).__init__(mu,
                                             num_loci,
                                             allele_len,
                                             rep,
                                             burnin,
                                             *args,
                                             **kwargs)

    def mutate(self, pop):
        dvars = pop.dvars()
        rep = dvars.rep
        gen = dvars.gen - self.burnin
        rng = sim.getRNG()
        mu = self.mu
        for i, ind in enumerate(pop.individuals()):
            for locus in range(self.num_loci):
                for ploidy in range(2):
                    if rng.randUniform() < mu[locus]:
                        try:
                            # At this step, idx holds intra-locus index
                            # of an available site.
                            idx = self.available[rep][locus].pop()
                        except IndexError:
                            if self.reclaim(pop, rep, locus):
                                idx = self.available[rep][locus].pop()
                            else:
                                # self.elongate(pop, rep, locus)
                                # idx = self.available[locus].pop()
                                sys.stderr.write(
                                    'rep={}, gen={}: '.format(rep, gen - burnin) +
                                    'maximum number of storable ' +
                                    'polymorphic sites reached')
                                return False
                        # We need to convert intra-locus index to
                        # inter-loci index starting from the beginning
                        # of a chromosome.
                        idx += locus * self.allele_len
                        ind.setAllele(1, idx, ploidy = ploidy)

        return True

    def reclaim(self, pop, rep, locus):
        '''Recycle sites, that are already fixed.'''
        allele_len = self.allele_len
        start = locus * allele_len
        stop = (locus + 1) * allele_len
        sites = range(start, stop)

        alleleStates = areAllelesMonomorphic(pop, sites)
        available = [i for i, state in enumerate(alleleStates) if state == True]

        if list(available) == list():
            return False

        for ind in pop.individuals():
            for site in available:
                ind.setAllele(0, site, ploidy = 0)
                ind.setAllele(0, site, ploidy = 1)
        self.available[rep][locus] = deque(available)
        return True

    def elongate(self, pop, locus):
        '''Increase the number of sites for a locus when no more site is avaialbe.'''
        raise NotImplementedError

# Implement the infinite-alleles model of mutation.  A storage for a
# locus is long int. When a new mutation arises, a new and unique
# value is assigned.  Therefore, we have to keep track of unused values.
class InfAlleleMutator(Mutator):

    def __init__(self, mu, num_loci, rep, *args, **kwargs):
        self.idx = [0] * num_loci
        super(InfAlleleMutator, self).__init__(mu,
                                               num_loci,
                                               1,
                                               rep,
                                               *args,
                                               **kwargs)

    def mutate(self, pop):
        dvars = pop.dvars()
        rep = dvars.rep
        gen = dvars.gen - self.burnin
        rng = sim.getRNG()
        mu = self.mu
        for i, ind in enumerate(pop.individuals()):
            for locus in range(self.num_loci):
                for ploidy in range(2):
                    if rng.randUniform() < mu[locus]:
                        self.idx[locus] += 1
                        ind.setAllele(self.idx[locus], locus, ploidy = ploidy)

        return True



def chunks(l, n):
    '''Divide a list into equal-length sublists'''
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def computeHeterozygosity(pop, num_loci, allele_len):
    # This only works if each locus has the same number of sites.
    h = [0] * num_loci                       # number of heterozygotes
    pop_size = float(pop.popSize())

    for ind in pop.individuals():
        genotype0 = chunks(ind.genotype(ploidy = 0), allele_len)
        genotype1 = chunks(ind.genotype(ploidy = 1), allele_len)

        for i, (g0, g1) in enumerate(zip(genotype0, genotype1)):
            if g0 != g1:
                h[i] += 1
    return list([i / pop_size for i in h])

def areAllelesMonomorphic(pop, alleleFreq = sim.ALL_AVAIL):
    sim.stat(pop, alleleFreq = alleleFreq)
    pop_size = pop.popSize()
    alleleNum = pop.dvars().alleleNum
    dip_size = 2 * pop_size
    return [False if alleleNum[allele][0] != 0 and alleleNum[allele][0] != dip_size else True
            for allele in alleleNum]


def computeNumberOfSegregatingSites(pop, num_loci, allele_len):
    alleleStates = areAllelesMonomorphic(pop)
    loci = chunks(alleleStates, allele_len)
    return list([len([True for site in locus if site == False]) for locus in loci])


# select unique a pair of parents, and make sure that father and
# mother are different individuals.
def pickTwoParents(pop):
    parents = list(pop.individuals())
    while True:
        yield random.sample(parents, 2)



if __name__ == '__main__':

    pop_size = args.POP
    ngen = args.NGEN
    nrep = args.NREP
    selfing_rate = args.selfing_rate
    recomb_rate = args.recombination_rate
    if recomb_rate == None:
        recomb_rate = 0.5
    else:
        recomb_rate /= (4. * pop_size)
        if recomb_rate > 0.5:
            sys.stderr.write('per-generation recombination rate > 0.5: use 0.5 instead\n')
            recomb_rate = 0.5
    num_loci = args.num_loci
    burnin = args.burn_in
    allele_len = args.num_segre_sites
    # convert scaled mutation rates to per-generation rates.
    mut_rates = [[pos, float(val) / (4. * pop_size)] for pos, val in args.mutation_rate]
    output = args.OUTPUT
    seed = args.seed
    if args.trajectory == None:
        to_explore = False
    else:
        to_explore = True
        traj_file = args.trajectory
    is_inf_alleles = args.infinite_alleles

    if seed > 0:
        sim.getRNG().set(seed = seed)

    # Define evolutionary operators here to avoid long lines later.

    if is_inf_alleles == False:
        population = sim.Population(size = pop_size,
                                    ploidy = 2,
                                    loci = num_loci * allele_len)
        mutator = InfSiteMutator(mu = mut_rates,
                                 num_loci = num_loci,
                                 allele_len = allele_len,
                                 rep = nrep,
                                 burnin = burnin)
        rec_sites = [i * allele_len for i in xrange(num_loci)]
        selfing = sim.SelfMating(ops = sim.Recombinator(rates = recomb_rate,
                                                        loci = rec_sites),
                                weight = pop_size * selfing_rate)
        offspring_func = sim.OffspringGenerator(
            ops = sim.Recombinator(rates = recomb_rate,
                                   loci = rec_sites))
        # write a header of result file here.  This is necessary as the
        # output is printed at each generation during exploration runs.
        if to_explore == True:
            writer = InfSiteWriter(traj_file, num_loci, allele_len, burnin)
            dump = [writer]
            writer.write_header()
        else:
            dump = []
    else:
        allele_len = 1
        population = sim.Population(size = pop_size,
                                    ploidy = 2,
                                    loci = num_loci * allele_len)
        mutator = InfAlleleMutator(mu = mut_rates,
                                   num_loci = num_loci,
                                   rep = nrep,
                                   burnin = burnin)
        selfing = sim.SelfMating(ops = sim.Recombinator(rates = recomb_rate),
                                 weight = pop_size * selfing_rate)
        offspring_func = sim.OffspringGenerator(ops = sim.Recombinator(rates = recomb_rate))
        # write a header of result file here.  This is necessary as the
        # output is printed at each generation during exploration runs.
        if to_explore == True:
            writer = InfAlleleWriter(traj_file, num_loci, burnin)
            dump = [writer]
            writer.write_header()
        else:
            dump = []

    parent_chooser = sim.PyParentsChooser(generator = pickTwoParents)
    outcross = sim.HomoMating(chooser = parent_chooser,
                              generator = offspring_func,
                              weight = pop_size * (1 - selfing_rate))

    # Population size is kept constant all the time.
    mating = sim.HeteroMating(matingSchemes = [selfing, outcross],
                              subPopSize = pop_size)

    simulator = sim.Simulator(pops = population,
                              rep = nrep)

    # Perform simulation.
    simulator.evolve(
        preOps = mutator,
        matingScheme = mating,
        postOps = dump,
        gen = ngen + burnin)


    # Because simuPOP only supports saving populations one-per-file, I
    # need to generate a bunch of file names.
    (path, ext) = os.path.splitext(output)
    if ext == '':                         # no extension
        filename = path + '.{}.pop'
    elif ext == '.pop':
        filename = path + '.{}' + ext
    else:
        filename = path + ext + '.{}.pop'

    # save the result if not in exploration runs.
    for pop in simulator.populations():
        if to_explore:
            writer.write(pop)
        pop.save(filename.format(pop.dvars().rep))
