# -*- mode: python; coding: utf-8; -*-

# infinite_sites.py - Implements simulations of partially linked pair
# of loci under the infinite sites model.

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

import simuOpt
simuOpt.setOptions(alleleType='binary')
import simuPop as simu
import partial_selfing.funcs.common as cf
import partial_selfing.funcs.infinite_sites as isf

def run(args):
    """Configure and run simulations."""
    pop = cf.get_population(size = orgs.NUM_IND,
                            loci = 2 * args.allele_length)


    mating_op = isf.get_mating_operator(r_rate = args.R_RATE,
                                        loci = 2,
                                        allele_length = args.allele_length,
                                        weight = args.S_RATE,
                                        size = args.NUM_IND)

    mating_op = isf.get_mutation_operator(m_rate = args.M_RATE,
                                          loci = 2,
                                          allele_length = args.allele_length,
                                          nrep = args.NUM_REP,
                                          burnin = args.burnin)

    simulator = simu.Simulator(pop = pop)

    simulator.evolve(
        preOps = mutation_op,
        matingScheme = mating_op,
        gen = args.NUM_GEN + args.burnin,
        rep = args.NUM_REP)
