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
        available = [start + i for i, state in enumerate(alleleStates) if state == True]

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


def computeHeterozygosity(pop, allele_len, num_loci):
    # This only works if each locus has the same number of sites.
    h = [0] * num_loci                       # number of heterozygotes
    pop_size = float(pop.popSize())

    for ind in pop.individuals():
        genotype0 = chunks(ind.genotype(ploidy = 0), allele_len)
        genotype1 = chunks(ind.genotype(ploidy = 1), allele_len)

        for i, (g0, g1) in enumerate(zip(genotype0, genotype1)):
            if g0 != g1:
                h[i] += 1
            else:
                print g0, g1
    return list([i / pop_size for i in h])

def areAllelesMonomorphic(pop, alleleFreq = sim.ALL_AVAIL):
    sim.state(pop, alleleFreq = alleleFreq)
    pop_size = pop.popSize()
    alleleNum = pop.dvars().alleleNum
    dip_size = 2 * pop_size
    return [False if alleleNum[allele][0] != 0 and alleleNum[allele][0] != dip_size else True
            for allele in alleleNum]


def computeNumberOfSegregatingSites(pop, allele_len, num_loci):
    alleleStates = areAllelesMonomorphic(pop)
    loci = chunks(alleleState, allele_len)
    return list([len([True for site in locus if site == False]) for locus in loci])


# select unique a pair of parents, and make sure that father and
# mother are different individuals.
def pickTwoParents(pop):
    parents = list(pop.individuals())
    while True:
        yield random.sample(parents, 2)


if __name__ == '__main__':

    if len(sys.argv) != 8:
        print('usage: inbreeding.py pop_size mut_rate0 mut_rate1 selfing_rate recomb_rate ngen nrep')
        sys.exit(1)

    pop_size = int(sys.argv[1])
    mut_rates = [float(d) for d in sys.argv[2:4]]
    selfing_rate, recomb_rate = [float(d) for d in sys.argv[4:6]]
    ngen, nrep = [int(i) for i in sys.argv[6:]]

    num_loci = 2                          # number of loci
    allele_len = 256                      # number of storable polymorphic sites per locus

    population = sim.Population(size = pop_size,
                                ploidy = 2,
                                loci = num_loci * allele_len)

    # Define evolutionary operators here to avoid long lines later.
    simulator = sim.Simulator(pops = population,
                              rep = nrep)

    genotype = sim.InitGenotype(genotype = [0, 0, 0, 0])

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


    # Perform simulation.
    simulator.evolve(
        initOps = genotype,
        preOps = mutator,
        matingScheme = mating,
        gen = ngen)

    # Extract heterozigosity (fraction of individuals that are
    # heterozigous) at the end of simulations.
    # Results are printed to STDOUT in CSV file.
    # print('"locus 0","locus 1"')
    for pop in simulator.populations():
        for ind in pop.individuals():
            print all([i == 0 for i in ind.genotype()])
        print computeHeterozygosity(pop, allele_len, num_loci)
        print computeNumberOfSegregatingSites(pop, allele_len, num_loci)
