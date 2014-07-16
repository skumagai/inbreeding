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


def get_init_genotype_by_prop(prop):
    s = sum(prop)
    return (len(prop), simu.InitGenotype(prop=[p / s for p in prop]))


def get_mutation_operator(m_rate, loci, nrep, burnin, new_idx=0):
    class MyMutator(simu.PyOperator):
        """
        A mutation operator class representing the infinite-alleles model.

        A new mutation is distinct from any other mutations already in a population.
        """
        def __init__(self):
            self.idx = list([new_idx] * loci for i in range(nrep))

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
                        if  rng.randUniform() < m_rate[locus]:
                            ind.setAllele(self.idx[rep][locus], locus, ploidy = ploidy)
                            self.idx[rep][locus] += 1
            return True

    return MyMutator()


def get_output_operator(config, field = 'self_gen'):
    output = config.outfile
    output_per = config.output_per
    N = config.N
    burnin = config.burnin
    ngen = config.gens

    header = [
        'replicate',
        'generation',
        'individual',
        'number of selfing',
        'chromosome'
    ] + ['locus {}'.format(i) for i in range(config.loci)])


    """Output genetic information of a population."""

    class MyWriter(simu.PyOperator):
        """A class handling output of genetic information of the entire population."""

        def __init__(self):
            with open(output, 'w') as f:
                writer = csv.DictWriter(f, header, delimiter = "\t")
                writer.writeheader()

            if output_per > 0:
                ats = [i + burnin for i in range(0, ngen, output_per)]
                super(MyWriter, self).__init__(func = self.write, at = ats)
            else:
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

                writer = csv.DictWriter(f, header, delimiter = "\t")

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
                                             [rep, gen, idx, selfing, ploidy] + geno)})

            return True

    return MyWriter()

def execute(config, pop, mating_op):

    init = config.init['model']

    if init == 'monomorphic':
        next_idx, init_genotype_op = cf.get_init_genotype_by_count(simu, 1)
    elif init == 'unique':
        next_idx, init_genotype_op = cf.get_init_genotype_by_count(simu, 2 * config.N)
    else:
        try:
            c = config.init['count']
            next_idx, init_genotype_op = cf.get_init_genotype_by_count(simu, c)
        except:
            next_idx, init_genotype_op = get_init_genotype_by_prop(config.init['freq'])

    init_info_op = cf.get_init_info(simu)

    mutation_op = get_mutation_operator(m_rate = config.m,
                                        loci = config.loci,
                                        nrep = config.reps,
                                        burnin = config.burnin,
                                        new_idx = next_idx)

    output_op = get_output_operator(config)

    simulator = simu.Simulator(pops = pop, rep = config.reps)

    if config.debug > 0:
        post_op = [simu.Stat(alleleFreq=simu.ALL_AVAIL, step=config.debug),
                   simu.PyEval(r"'%s\n' % alleleFreq", step=config.debug)]
    else:
        post_op = []

    if config.output_per > 0:
        post_op.append(output_op)

    simulator.evolve(
        initOps = [init_info_op, init_genotype_op],
        preOps = mutation_op,
        matingScheme = mating_op,
        postOps = post_op,
        finalOps = output_op,
        gen = config.gens + config.burnin)



def run(config):
    if config.model == 'androdioecy':
        cf.androdioecy(simu, execute, config)
    elif config.model == 'gynodioecy':
        cf.gynodioecy(simu, execute, config)
    else:
        cf.pure_hermaphrodite(simu, execute, config)
