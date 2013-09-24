# -*- mode: python; coding: utf-8; -*-

# infinite_loci.py - Implements simulations of partially linked pair
# of loci under the infinite alleles model.

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

import csv

import simuOpt
simuOpt.setOptions(alleleType='long')
import simuPOP as simu
import partial_selfing.common as cf


def get_mating_operator(r_rate, weight, size, field='self_gen'):
    """
    Construct partially selfing mating operator under the infinite alleles model.

    A fraction, 0 <= weight <= 1, of offspring is generated by selfing, and others are
    generated by outcrossing.  In this model, there is no specific sex so that any
    individual can mate with any other individuals in a population.
    Furthermore, a parent can participate in both selfing and outcrossing.
    """
    selfing = simu.SelfMating(ops = [simu.Recombinator(rates = r_rate),
                                     cf.MySelfingTagger(field)],
                              weight = weight)

    outcross = simu.HomoMating(chooser = simu.PyParentsChooser(generator = cf.pickTwoParents),
                               generator = simu.OffspringGenerator(
                                   ops = [simu.Recombinator(rates = r_rate),
                                          cf.MyOutcrossingTagger(field)]),
                               weight = 1.0 - weight)
    return simu.HeteroMating(matingSchemes = [selfing, outcross],
                             subPopSize = size)


def get_mutation_operator(m_rate, loci, nrep, burnin):
    class MyMutator(simu.PyOperator):
        """
        A mutation operator class representing the infinite-alleles model.

        A new mutation is distinct from any other mutations already in a population.
        """
        def __init__(self):
            self.idx = list([0] * loci for i in range(nrep))

            super(MyMutator, self).__init__(func = self.mutate)


        def mutate(self, pop):
            """Add mutations to organisms."""
            rng = simu.getRNG()

            dvars = pop.dvars()
            rep = dvars.rep
            gen = dvars.gen - burnin

            for i, ind in enumerate(pop.individuals()):
                for locus in range(loci):
                    for ploidy in range(2):
                        if rng.randUniform() < m_rate[locus]:
                            self.idx[rep][locus] += 1
                            ind.setAllele(self.idx[rep][locus], locus, ploidy = ploidy)
            return True

    return MyMutator()


def get_output_operator(size,
                        ngen,
                        nrep,
                        m_rate,
                        r_rate,
                        s_rate,
                        loci,
                        burnin,
                        output,
                        field = 'self_gen'):

    data = ['infinite alleles',
            size,
            ngen,
            nrep,
            loci,
            m_rate,
            s_rate,
            r_rate,
            burnin,
            simu.getRNG().seed()]

    header = ['mutation model',
              'number of individuals',
              'number of generations',
              'number of replicates',
              'number of loci',
              'mutation rate',
              'selfing rate',
              'recombination rate',
              'number of burnin generations',
              'random number seed',
              'replicate',
              'generation',
              'individual',
              'number of selfing',
              'chromosome'] + ['locus {}'.format(i) for i in range(loci)]


    """Output genetic information of a population."""

    class MyWriter(simu.PyOperator):
        """A class handling output of genetic information of the entire population."""

        def __init__(self):
            with open(output, 'w') as f:
                writer = csv.DictWriter(f, header)
                writer.writeheader()

            super(MyWriter, self).__init__(func = self.write)


        def write(self, pop):
            # In order to keep output file structure simple, all
            # information regarding to simulation such as model
            # parameters are included into each row.  This obviously
            # caused repetition of simulation-wide parameters many
            # times and excessive use of storage space.  However, I
            # consider an upside, the simplicity of the output file
            # structure, is well worth the cost.

            with open(output, 'a') as f:
                dvars = pop.dvars()
                rep = dvars.rep
                gen = dvars.gen - burnin

                writer = csv.DictWriter(f, header)

                # write genotype row by row.  Each row contains a list
                # of genes on a single chromosome.  Because simulated
                # organisms are diploid, each individual occupy two
                # (successive) rows.
                for idx, ind in enumerate(pop.individuals()):
                    selfing = ind.info(field)
                    for ploidy in range(2):
                        geno = list(ind.genotype(ploidy = ploidy))
                        writer.writerow({key: value for key, value in
                                         zip(header,
                                             data + [rep, gen, idx, selfing, ploidy] + geno)})

            return True

    return MyWriter()


def run(args):
    pop = cf.get_population(size=args.NUM_IND,
                            loci = args.NUM_LOCI)

    init_info_op = cf.get_init_info()

    init_genotype_op = cf.get_init_genotype()

    mating_op = get_mating_operator(r_rate=args.R_RATE,
                                    weight = args.S_RATE,
                                    size = args.NUM_IND)

    mutation_op = get_mutation_operator(m_rate = args.M_RATE,
                                        loci = args.NUM_LOCI,
                                        nrep = args.NUM_REP,
                                        burnin = args.burnin)

    output_op = get_output_operator(size = args.NUM_IND,
                                    m_rate = args.M_RATE,
                                    r_rate = args.R_RATE,
                                    s_rate = args.S_RATE,
                                    loci = args.NUM_LOCI,
                                    nrep = args.NUM_REP,
                                    ngen = args.NUM_GEN,
                                    burnin = args.burnin,
                                    output = args.OUTFILE)


    simulator = simu.Simulator(pops = pop, rep = args.NUM_REP)


    simulator.evolve(
        initOps = [init_info_op, init_genotype_op],
        preOps = mutation_op,
        matingScheme = mating_op,
        finalOps = output_op,
        gen = args.NUM_GEN + args.burnin)
