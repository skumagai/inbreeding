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

import sys
import random

import simuOpt
simuOpt.setOptions(alleleType = 'long',
                   optimized = False,
                   quiet = True)

import simuPOP as sim

# Implement infinite-sites mutation model.  Whenever a new mutation
# arises, assign n + 1-th allele to a mutatnt, for n the maximum
# observed allele.  Alleles are encoded as non-negative integer, where
# zero indicates ancestral type.  The largest allele indicates the
# number of mutations occurred to each locus
class InfSiteMutator(sim.PyOperator):

    def __init__(self, mu, pop_size, *args, **kwargs):
        self.mu = mu
        self.idx = [0, 0]
        self.pop_size = pop_size
        super(InfSiteMutator, self).__init__(func = self.mutate, *args, **kwargs)

    def mutate(self, pop):
        rng = sim.getRNG()
        (mu0, mu1) = self.mu
        (idx0, idx1) = self.idx
        rngs = [rng.randUniform() for i in range(4 * self.pop_size)]
        for i, ind in enumerate(pop.individuals()):
            # manually unrolling inner loops to (hopefully) speed up
            # simulation.

            # The first two elements are for the first locus (both chromosomes)
            if rngs[i] < mu0:
                idx0 += 1
                ind.setAllele(idx0, 0)

            if rngs[i+1] < mu0:
                idx0 += 1
                ind.setAllele(idx0, 2)

            # The last two elements are for the second locus (both chromosomes)
            if rngs[i+2] < mu1:
                idx1 += 1
                ind.setAllele(idx1, 1)

            if rngs[i+3] < mu1:
                idx1 += 1
                ind.setAllele(idx1, 3)

        self.idx = [idx0, idx1]
        return True

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

    population = sim.Population(size = pop_size,
                                ploidy = 2,
                                loci = 2)

    # Define evolutionary operators here to avoid long lines later.
    simulator = sim.Simulator(pops = population,
                              rep = nrep)

    genotype = sim.InitGenotype(genotype = [0, 0, 0, 0])

    mutator = InfSiteMutator(mut_rates, pop_size)

    # (selfing_rate * 100) % of individuals are selfers.
    selfing = sim.SelfMating(ops = sim.Recombinator(rates = recomb_rate),
                             weight = pop_size * selfing_rate)

    # the rest are outcrossers.
    offspring_func = sim.OffspringGenerator(ops = sim.Recombinator(rates = recomb_rate))

    parent_chooser = sim.PyParentsChooser(generator = pickTwoParents)
    outcross = sim.HomoMating(chooser = parent_chooser,
                              generator = offspring_func,
                              weight = pop_sizes * (1 - selfing_rate))

    # Population size is kept constant all the time.
    mating = sim.HeteroMating(matingSchemes = [selfing, outcross],
                              subPopSize = pop_size)


    # Perform simulation.
    simulator.evolve(
        initOps = genotype,
        preOps = mutator,
        matingScheme = mating,
        finalOps = sim.Stat(heteroFreq = [0,1]),
        gen = ngen)

    # Extract heterozigosity (fraction of individuals that are
    # heterozigous) at the end of simulations.
    # Results are printed to STDOUT in CSV file.
    print('"locus 0","locus 1"\n')
    for pop in simulator.populations():
        vars = pop.dvars()
        print('{},{}\n'.format(vars.heteroFreq[0], vars.heteroFreq[1])),
