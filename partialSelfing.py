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

import simuOpt
simuOpt.setOptions(alleleType = 'binary',
                   optimized = False,
                   quiet = True)

import simuPOP as sim

# Implement infinite-sites mutation model.  As new mutations arise,
# number of polymorphic sites increases.  In the current
# implementation a fix number of slots are pre-allocated to store
# polymorphic sites.  Once all available space is exhausted, sites
# that are monomorphic (due to random drift) are recycled.  When the
# recycle is failed to open up any new space, simulations are
# terminated.

class Store(sim.PyOperator):

    def __init__(self, num_loci, allele_len, rep, gen, *args, **kwargs):
        self.num_loci = num_loci
        self.allele_len = allele_len
        self.h = {r: {n: [0] * num_loci for n in range(gen)}
                  for r in range(rep)}
        self.segre = [[[0] * num_loci] * gen] * rep
        super(Store, self).__init__(func = self.store, *args, **kwargs)

    def store(self, pop):
        dvars = pop.dvars()
        rep = dvars.rep
        gen = dvars.gen
        num_loci = self.num_loci
        allele_len = self.allele_len
        for i, val in enumerate(computeHeterozygosity(pop, num_loci, allele_len)):
            self.h[rep][gen][i] = val
        for i, val in enumerate(computeNumberOfSegregatingSites(pop, num_loci, allele_len)):
            self.segre[rep][gen][i] = val
        return True



#
# Dynamically increasing the number of slots is planned but not
# implemented.
class InfSiteMutator(sim.PyOperator):

    def __init__(self, mu, num_loci, allele_len, rep, *args, **kwargs):
        self.mu = mu
        self.num_loci = num_loci
        self.allele_len = allele_len
        self.available = {r: [deque(range(allele_len)) for i in range(2)] for r in range(rep)}
        super(InfSiteMutator, self).__init__(func = self.mutate, *args, **kwargs)

    def mutate(self, pop):
        rep = pop.dvars().rep
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
                                sys.exit('maximum number of storable polymorphic sites reached')
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
                        nargs='?',
                        type=argparse.FileType('w'),
                        default=sys.stdout,
                        help='output file name (default: STDOUT)')
    parser.add_argument('-r', '--recombination-rate',
                        metavar='R',
                        type=float,
                        default=0.5,
                        help='recombination rate (default: 0.5)')
    parser.add_argument('-s', '--selfing-rate',
                        metavar='S',
                        type=float,
                        default=0,
                        help='selfing rate (default: 0)')
    parser.add_argument('-m', '--mutation-rate',
                        metavar='M',
                        type=float,
                        nargs='*',
                        default=[0, 0],
                        help='mutation rate (default:0)')
    parser.add_argument('-b', '--burn-in',
                        metavar='B',
                        type=int,
                        default=0,
                        help='burn-in (default: 0)')
    parser.add_argument('-n', '--num_loci',
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
    parser.add_argument('--explore',
                        action='store_true',
                        help='record heterozygosity and number of segregating sites each generation for later inspection to determine an appropriate durtion of burn-in')
    return parser.parse_args()



if __name__ == '__main__':

    # if len(sys.argv) != 8:
    #     print('usage: inbreeding.py pop_size mut_rate0 mut_rate1 selfing_rate recomb_rate ngen nrep')
    #     sys.exit(1)

    args = parseArgs()
    pop_size = args.POP
    ngen = args.NGEN
    nrep = args.NREP
    selfing_rate = args.selfing_rate
    recomb_rate = args.recombination_rate
    num_loci = args.num_loci
    burnin = args.burn_in
    allele_len = args.num_segre_sites
    mut_rates = args.mutation_rate
    output = args.OUTPUT
    seed = args.seed
    to_explore = args.explore

    if seed > 0:
        sim.getRNG().set(seed = seed)


    if len(mut_rates) != 1 and len(mut_rates) != num_loci:
        sys.exit('number of mutation rates must be 1 or equal to the number of loci')

    population = sim.Population(size = pop_size,
                                ploidy = 2,
                                loci = num_loci * allele_len)

    # Define evolutionary operators here to avoid long lines later.
    simulator = sim.Simulator(pops = population,
                              rep = nrep)

    mutator = InfSiteMutator(mu = mut_rates,
                             num_loci = num_loci,
                             allele_len = allele_len,
                             rep = nrep)

    # (selfing_rate * 100) % of individuals are selfers.
    selfing = sim.SelfMating(ops = sim.Recombinator(rates = recomb_rate),
                             weight = pop_size * selfing_rate)

    # the rest are outcrossers.
    offspring_func = sim.OffspringGenerator(ops = sim.Recombinator(rates = recomb_rate))

    parent_chooser = sim.PyParentsChooser(generator = pickTwoParents)
    outcross = sim.HomoMating(chooser = parent_chooser,
                              generator = offspring_func,
                              weight = pop_size * (1 - selfing_rate))

    # Population size is kept constant all the time.
    mating = sim.HeteroMating(matingSchemes = [selfing, outcross],
                              subPopSize = pop_size)

    store = Store(num_loci, allele_len, nrep, ngen)

    # Perform simulation.
    simulator.evolve(
        preOps = mutator,
        matingScheme = mating,
        postOps = store,
        gen = ngen)

    # Extract heterozigosity (fraction of individuals that are
    # heterozigous) at the end of simulations.
    # Results are printed to STDOUT in CSV file.
    # print('"locus 0","locus 1"')
    output.write('"rep","gen","1-f (locus 0)","1-f (locus1)","# segre (locus 0)", "# segre (locus 1)"\n')
    if to_explore == True:
        for r in xrange(nrep):
            for n in xrange(ngen):
                output.write('{},{},{},{},{},{}\n'.format(r,
                                                          n,
                                                          *(store.h[r][n] + store.segre[r][n])))
    else:
        for pop in simulator.populations():
            dvars = pop.dvars()
            vals = computeHeterozygosity(pop, num_loci, allele_len) + \
                   computeNumberOfSegregatingSites(pop, num_loci, allele_len)
            output.write("{},{},{},{},{},{}\n".format(dvars.rep,
                                                      dvars.gen,
                                                      *vals))
