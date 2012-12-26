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
import json
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
    parser.add_argument('-f', '--force-replication',
                        action='store_true',
                        help='instead of using several chromosomes to represents replication, run several separate simulations.  Note that this can be much slower than using chromosomes.')
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

    def __init__(self, output, num_loci, allele_len, rep_mode, burnin, *args, **kwargs):
        self.num_loci = num_loci
        self.allele_len = allele_len
        self.output = output
        self.burnin = burnin
        self.has_header_printed = False
        self.rep_mode = rep_mode
        super(Writer, self).__init__(func = self.write, *args, **kwargs)

    def write(self, pop):
        '''stub for a function to write data'''
        raise NotImplementedError

    def write_header(self):
        '''stub for a function to write header'''
        raise NotImplementedError

    def _write_common(self, pop, rep, gen, chrom = 0):
        self.output.write('{},{},'.format(rep, gen - self.burnin))
        self.output.write(','.join(['{}'.format(val)
                                    for val in computeHeterozygosity(pop,
                                                                     self.num_loci,
                                                                     self.allele_len,
                                                                     chrom)]))

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

    def __init__(self, output, num_loci, allele_len, rep_mode, burnin, *args, **kwargs):
        if rep_mode is True:
            self.write = self.writeRep
        else:
            self.write = self.writeChrom

        super(InfSiteWriter, self).__init__(output,
                                            num_loci,
                                            allele_len,
                                            rep_mode,
                                            burnin,
                                            *args,
                                            **kwargs)

    def writeRep(self, pop):
        dvars = pop.dvars()
        self._write_common(pop, dvars.rep, dvars.gen)
        self.output.write(','.join(['{}' for i in xrange(self.num_loci)]).format(
            *computeNumberOfSegregatingSites(pop, num_loci, allele_len)) + '\n')
        return True

    def writeChrom(self, pop, rep):
        gen = pop.dvars().gen
        for rep in xrange(pop.numChrom()):
            self._write_common(pop, rep, gen, rep)
            self.output.write(','.join(['{}' for i in xrange(self.num_loci)]).format(
                *computeNumberOfSegregatingSites(pop, num_loci, allele_len, rep)) + '\n')
        return True


    def write_header(self):
        self._write_common_header()
        self.output.write(',' + ','.join(['"# segre (locus {})"'.format(i)
                                          for i in xrange(self.num_loci)]) + '\n')


class InfAlleleWriter(Writer):
    '''Write heterozigosities at each generation.

    This class is intended to be used in explorative runs and only
    under the infinite-alleles model.'''

    def __init__(self, output, num_loci, allele_len, rep_mode, burnin, *args, **kwargs):
        if rep_mode is True:
            self.write = self.writeRep
        else:
            self.write = self.writeChrom

        super(InfAlleleWriter, self).__init__(output,
                                              num_loci,
                                              1,
                                              rep_mode,
                                              burnin,
                                              *args,
                                              **kwargs)

    def writeRep(self, pop):
        self._write_common(pop)
        self.output.write('\n')
        return True

    def writeChrom(self, pop):
        for rep in xrange(pop.numChrom()):
            self._write_common(pop, rep)
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
    def __init__(self, mu, num_loci, allele_len, rep, rep_mode, burnin, *args, **kwargs):
        self.num_loci = num_loci
        self._build_mutation_rates(mu)
        self.allele_len = allele_len
        self.burnin = burnin
        self.rep = rep
        self.rep_mode = rep_mode
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

    def __init__(self, mu, num_loci, allele_len, rep, rep_mode, burnin, *args, **kwargs):
        self.available = {r: [deque(xrange(allele_len))
                              for i in range(num_loci)]
                          for r in range(rep)}
        if rep_mode is True:
            self.mutate = self.mutateRep
        else:
            self.mutate = self.mutateChrom

        super(InfSiteMutator, self).__init__(mu,
                                             num_loci,
                                             allele_len,
                                             rep,
                                             rep_mode,
                                             burnin,
                                             *args,
                                             **kwargs)

    def mutateRep(self, pop):
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


    def mutateChrom(self, pop):
        dvars = pop.dvars()
        gen = dvars.gen - self.burnin
        rng = sim.getRNG()
        mu = self.mu
        for i, ind in enumerate(pop.individuals()):
            for locus in xrange(self.num_loci):
                for ploidy in xrange(2):
                    for rep in xrange(self.rep):
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
                            ind.setAllele(1, idx, ploidy = ploidy, chrom = rep)

        return True


    # TODO: adjust to multiple chromosomes (for replications)
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

    def __init__(self, mu, num_loci, allele_len, rep, rep_mode, *args, **kwargs):
        self.idx = [[0] * num_loci for i in xrange(rep)]
        if rep_mode is True:
            self.mutate = self.mutateRep
        else:
            self.mutate = self.mutateChrom

        super(InfAlleleMutator, self).__init__(mu,
                                               num_loci,
                                               1,
                                               rep,
                                               rep_mode
                                               *args,
                                               **kwargs)

    def mutateRep(self, pop):
        dvars = pop.dvars()
        rep = dvars.rep
        gen = dvars.gen - self.burnin
        rng = sim.getRNG()
        mu = self.mu
        for i, ind in enumerate(pop.individuals()):
            for locus in range(self.num_loci):
                for ploidy in range(2):
                    if rng.randUniform() < mu[locus]:
                        self.idx[rep][locus] += 1
                        ind.setAllele(self.idx[rep][locus], locus, ploidy = ploidy)

        return True


    def mutateChrom(self, pop):
        dvars = pop.dvars()
        gen = dvars.gen - self.burnin
        rng = sim.getRNG()
        mu = self.mu
        for i, ind in enumerate(pop.individuals()):
            for locus in xrange(self.num_loci):
                for ploidy in xrange(2):
                    for rep in xrange(self.rep):
                        if rng.randUniform() < mu[locus]:
                            self.idx[rep][locus] += 1
                            ind.setAllele(self.idx[rep][locus],
                                          locus,
                                          ploidy = ploidy,
                                          chrom = rep)

        return True



def chunks(l, n):
    '''Divide a list into equal-length sublists'''
    for i in xrange(0, len(l), n):
        yield l[i:i+n]


def computeHeterozygosity(pop, num_loci, allele_len, chrom):
    # This only works if each locus has the same number of sites.
    h = [0] * num_loci                       # number of heterozygotes
    pop_size = float(pop.popSize())

    for ind in pop.individuals():
        genotype0 = chunks(ind.genotype(ploidy = 0, chrom = chrom), allele_len)
        genotype1 = chunks(ind.genotype(ploidy = 1, chrom = chrom), allele_len)

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


def computeNumberOfSegregatingSites(pop, num_loci, allele_len, chrom = 0):
    alleleStates = areAllelesMonomorphic(pop)
    start = chrom * num_loci * allele_len
    stop = start + num_loci * allele_len
    loci = chunks(alleleStates[start:stop], allele_len)
    return list([len([True for site in locus if site == False]) for locus in loci])


# select unique a pair of parents, and make sure that father and
# mother are different individuals.
def pickTwoParents(pop):
    parents = list(pop.individuals())
    while True:
        yield random.sample(parents, 2)

def scaleParam(param, pop_size):
    return 4. * pop_size * float(param)

def unscaleParam(param, pop_size):
    return float(param) / (4. * pop_size)

def buildOutputFileName(output):
    (path, ext) = os.path.splitext(output)
    subst = '.{:0{width}d}'
    if ext == '':                         # no extension
        filename = path + subst + '.pop'
    elif ext == '.pop':
        filename = path + subst + ext
    else:
        filename = path + ext + subst + '.pop'

    return filename

def main():

    pop_size = args.POP
    ngen = args.NGEN
    nrep = args.NREP
    selfing_rate = args.selfing_rate
    recomb_rate = args.recombination_rate
    if recomb_rate == None:
        recomb_rate = 0.5
    else:
        recomb_rate = unscaleParam(recomb_rate, pop_size)
        if recomb_rate > 0.5:
            sys.stderr.write('per-generation recombination rate > 0.5: use 0.5 instead\n')
            recomb_rate = 0.5
    num_loci = args.num_loci
    burnin = args.burn_in
    allele_len = args.num_segre_sites

    # convert scaled mutation rates to per-generation rates.
    mut_rates = [[pos, unscaleParam(val, pop_size)] for pos, val in args.mutation_rate]
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

    if is_inf_alleles == True:
        allele_len = 1
        mutator = InfAlleleMutator
        writer = InfAlleleWriter
        rec_sites = sim.ALL_AVAIL
    else:
        mutator = InfSiteMutator
        writer = InfSiteWriter
        rec_sites = [i * allele_len for i in xrange(num_loci)]

    # construct a blue-print of a population.
    if args.force_replication is True:
        population = sim.Population(size = pop_size,
                                    ploidy = 2,
                                    loci = num_loci * allele_len)
    else:
        population = sim.Population(size = pop_size,
                                    ploidy = 2,
                                    loci = [num_loci * allele_len] * nrep)
        # run only a single run of a simulation
        nrep = 1


    # Define operators used during simulations.
    mutator = mutator(mu = mut_rates,
                      num_loci = num_loci,
                      allele_len = allele_len,
                      rep = nrep,
                      rep_mode = args.force_replication,
                      burnin = burnin)
    selfing = sim.SelfMating(ops = sim.Recombinator(rates = recomb_rate,
                                                    loci = rec_sites),
                                                    weight = pop_size * selfing_rate)
    offspring_func = sim.OffspringGenerator(
        ops = sim.Recombinator(rates = recomb_rate,
                               loci = rec_sites))

    parent_chooser = sim.PyParentsChooser(generator = pickTwoParents)
    outcross = sim.HomoMating(chooser = parent_chooser,
                              generator = offspring_func,
                              weight = pop_size * (1 - selfing_rate))

    # Population size is kept constant all the time.
    mating = sim.HeteroMating(matingSchemes = [selfing, outcross],
                              subPopSize = pop_size)

    # Because simuPOP only supports saving populations one-per-file, I
    # need to generate a bunch of file names.
    filename = buildOutputFileName(output)
    width = str(len(str(nrep - 1)))
    finalops = [sim.SavePopulation('!"' + filename + '".format(rep, width=' + width + ')')]

    if to_explore == True:
        traj_writer = writer(traj_file, num_loci, allele_len, burnin)
        finalops.append(traj_writer)
        # write a header of result file here.  This is necessary as the
        # output is printed at each generation during exploration runs.
        traj_writer.write_header()
        dump = traj_writer
    else:
        dump = []


    simulator = sim.Simulator(pops = population,
                              rep = nrep)

    # Perform simulation.
    simulator.evolve(
        preOps = mutator,
        matingScheme = mating,
        postOps = dump,
        finalOps = finalops,
        gen = ngen + burnin)


    # save information about simulations to a file.
    with open('conf.json', 'w') as wf:
        if args.infinite_alleles:
            mmode = u'infinite-alleles'
        else:
            mmode = u'infinite-sites'
        files = [filename.format(i, width=width) for i in xrange(nrep)]

        if args.force_replication:
            rmode = u'by chromosomes'
        else:
            rmode = u'by repeated simulations'

        if len(set(mutator.mu)) == 1:
            m = mutator.mu[0]
            sm = scaleParam(m, pop_size)
        else:
            m = mutator.mu
            sm = [scaleParam(mm, pop_size) for mm in m]

        info = {u'mutation rate': {u'unscaled': m,
                                   u'scaled': sm},
                u'recombination rate': {u'unscaled': recomb_rate,
                                        u'scaled': scaleParam(recomb_rate, pop_size)},
                u'selfing rate': selfing_rate,
                u'population size': pop_size,
                u'number of replications': nrep,
                u'replication mode': rmode,
                u'mutation mode': mmode,
                u'number of loci': num_loci,
                u'files': files,
                u'seed': hex(sim.getRNG().seed())}

        json.dump(info, wf)

if __name__ == '__main__':
    main()
