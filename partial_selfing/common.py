# -*- mode: python; coding: utf-8; -*-

# common.py - common functions usable under different mutational models.

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

import simuPOP.sampling as sampling

def get_population(simu, size, loci, infoFields='self_gen'):
    """Construct a population object."""
    return simu.Population(size = size,
                           ploidy = 2,
                           loci = loci,
                           infoFields = infoFields)


def get_init_info(simu, field='self_gen'):
    """Zero initialize info field `field`."""
    return simu.InitInfo(0, infoFields=field)


def get_init_genotype_by_count(simu, n):
    """
    Set genotype of inital population.
    """
    return (n, simu.InitGenotype(prop=[1.0 / n for dummy in range(n)]))


def pick_pure_hermaphrodite_parents(simu, a, tau):
    rng = simu.getRNG()
    runif = rng.randUniform
    rint = rng.randInt
    def generator(pop):
        N = pop.popSize()
        while True:
            if runif() < a:         # uniparental
                if runif() < tau:   # zygote survived
                    # print(1)
                    yield pop.individual(rint(N))
            else:                   # biparental
                pair = [rint(N), rint(N)]
                while pair[0] == pair[1]:
                    pair[1] = rint(N)
                # print(2)
                yield [pop.individual(p) for p in pair]
    return generator


def pick_androdioecious_parents(simu, a, tau, sigma):
    rng = simu.getRNG()
    runif = rng.randUniform
    rint = rng.randInt
    def generator(pop):
        gen = -1
        while True:
            ngen = pop.dvars().gen
            if gen != ngen:
                # At the beginning of a generation, extract the
                # sex-specific subpopulations from a parental
                # population. The sex-specific subpopulations are used
                # throughout mating events in one generation.
                gen = ngen
                m = pop.extractSubPops(subPops = [(0, 0)])
                h = pop.extractSubPops(subPops = [(0, 1)])
                Nm = m.popSize()
                Nh = h.popSize()

            if runif() < a:         # uniparental
                if runif() < tau:   # zygote survives
                    # print(1)
                    yield h.individual(rint(Nh))
            else:                   # biparental
                if runif() < sigma: # successful fertilization
                    # print(2)
                    yield [m.individual(rint(Nm)), h.individual(rint(Nh))]
    return generator


def pick_gynodioecious_parents(simu, a, tau, sigma):
    rng = simu.getRNG()
    runif = rng.randUniform
    rint = rng.randInt
    def generator(pop):
        gen = -1
        while True:
            ngen = pop.dvars().gen
            if gen != ngen:
                # At the beginning of a generation, extract the
                # sex-specific subpopulations from a parental
                # population. The sex-specific subpopulations are used
                # throughout mating events in one generation.
                gen = ngen
                h = pop.extractSubPops(subPops = [(0, 0)])
                f = pop.extractSubPops(subPops = [(0, 1)])
                Nh = h.popSize()
                Nf = f.popSize()

            if runif() < Nh / (Nh + Nf * sigma): # the seed is from hermaphrodite
                if runif() < a:                  # uniparental
                    if runif() < tau:            # zygote survived
                        # print(1)
                        yield h.individual(rint(Nh))
                else:           # biparental
                    first, second = rint(Nh), rint(Nh)
                    while first == second:
                        second = rint(Nh)
                    # print(3)
                    yield [h.individual(first), h.individual(second)]
            else:               # the seed is from female.
                # print(2)
                yield [h.individual(rint(Nh)), f.individual(rint(Nf))]
    return generator


def get_selfing_tagger(simu, field):
    class MySelfingTagger(simu.PyOperator):
        """
        Update information field to reflect selfing.

        When selfing occurred, this operator record the fact by incrementing the value
        of `field` by one.
        """

        def __init__(self, field='self_gen'):
            self.field = field
            super(MySelfingTagger, self).__init__(func = self.record)

        def record(self, pop, off, dad, mom):
            """
            Increment the value of offspring's infofield `field` by one if it is uniparental.
            Otherwise reset the value to 0.
            """
            if mom is not None:
                off.setInfo(0, self.field)
            else:
                off.setInfo(dad.info(self.field) + 1, self.field)
            return True
    return MySelfingTagger(field)

def get_pure_hermaphrodite_mating(simu, r_rate, a, tau, size, rec_sites, field='self_gen'):
    """
    Construct mating scheme for pure hermaphrodite with partial selfing under
    the infinite alleles model.

    A fraction, 0 <= weight <= 1, of offspring is generated by selfing, and others are
    generated by outcrossing.  In this model, there is no specific sex so that any
    individual can mate with any other individuals in a population.
    Furthermore, a parent can participate in both selfing and outcrossing.
    """

    parents_chooser = simu.PyParentsChooser(
        pick_pure_hermaphrodite_parents(simu = simu, a = a, tau = tau)
    )

    selfing_tagger = get_selfing_tagger(simu, field)

    return simu.HomoMating(chooser = parents_chooser,
                           generator = simu.OffspringGenerator(
                               ops = [simu.Recombinator(rates = r_rate, loci = rec_sites),
                                      selfing_tagger]),
                           subPopSize = size)


def get_androdioecious_mating(simu, r_rate, a, tau, sigma,
                              size, sex_ratio, rec_sites, field = 'self_gen'):

    sexMode = (simu.PROB_OF_MALES, sex_ratio)

    parents_chooser = simu.PyParentsChooser(
        pick_androdioecious_parents(simu = simu, a = a, tau = tau, sigma = sigma)
    )

    selfing_tagger = get_selfing_tagger(simu, field)
    return simu.HomoMating(chooser = parents_chooser,
                           generator = simu.OffspringGenerator(
                               ops = [simu.Recombinator(rates = r_rate, loci = rec_sites),
                                      selfing_tagger],
                               sexMode = sexMode),
                           subPopSize = size)


def get_gynodioecious_mating(simu, r_rate, a, tau, sigma,
                             size, sex_ratio, rec_sites, field = 'self_gen'):

    sexMode = (simu.PROB_OF_MALES, sex_ratio)

    parents_chooser = simu.PyParentsChooser(
        pick_gynodioecious_parents(simu = simu, a = a, tau = tau, sigma = sigma)
    )

    selfing_tagger = get_selfing_tagger(simu, field)
    return simu.HomoMating(chooser = parents_chooser,
                           generator = simu.OffspringGenerator(
                               ops = [simu.Recombinator(rates = r_rate, loci = rec_sites),
                                      selfing_tagger],
                               sexMode = sexMode),
                           subPopSize = size)


def pure_hermaphrodite(simu, execute_func, config):
    pop = get_population(simu = simu,
                         size = config.N,
                         loci = config.loci * config.allele_length)

    # Index of sites, after which recombinations happen.
    rec_loci = [config.allele_length * i - 1 for i in range(1, config.loci)]

    mating_op = get_pure_hermaphrodite_mating(simu,
                                              r_rate = config.r,
                                              a = config.a,
                                              tau = config.tau,
                                              size = config.N,
                                              rec_sites = rec_loci)

    execute_func(config, pop, mating_op)


def androdioecy(simu, execute_func, config):
    pop = get_population(simu = simu,
                         size = config.N,
                         loci = config.loci * config.allele_length)

    sex_ratio = 1.0 - config.sex_ratio
    simu.initSex(pop, maleFreq = sex_ratio)
    pop.setVirtualSplitter(simu.SexSplitter())

    # Index of sites, after which recombinations happen.
    rec_loci = [config.allele_length * i - 1 for i in range(1, config.loci)]

    mating_op = get_androdioecious_mating(simu,
                                          r_rate = config.r,
                                          a = config.a,
                                          tau = config.tau,
                                          sigma = config.sigma,
                                          size = config.N,
                                          sex_ratio = sex_ratio,
                                          rec_sites = rec_loci)

    execute_func(config, pop, mating_op)



def gynodioecy(simu, execute_func, config):
    pop = get_population(simu = simu,
                         size = config.N,
                         loci = config.loci * config.allele_length)

    sex_ratio = config.sex_ratio
    simu.initSex(pop, maleFreq = sex_ratio)
    pop.setVirtualSplitter(simu.SexSplitter())

    # Index of sites, after which recombinations happen.
    rec_loci = [config.allele_length * i - 1 for i in range(1, config.loci)]

    mating_op = get_gynodioecious_mating(simu,
                                         r_rate = config.r,
                                         a = config.a,
                                         tau = config.tau,
                                         sigma = config.sigma,
                                         size = config.N,
                                         sex_ratio = sex_ratio,
                                         rec_sites = rec_loci)

    execute_func(config, pop, mating_op)
