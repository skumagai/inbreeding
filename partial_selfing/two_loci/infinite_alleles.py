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

import simuOpt
simuOpt.setOptions(alleleType='long')
import simuPOP as simu
import partial_selfing.funcs.common as cf
import partial_selfing.funcs.infinite_alleles as iaf

def run(args):
    pop = cf.get_population(size=args.NUM_IND,
                            loci = 2)

    init_info_op = cf.get_init_info()

    init_genotype_op = cf.get_init_genotype()

    mating_op = iaf.get_mating_operator(r_rate=args.R_RATE,
                                        weight = args.S_RATE,
                                        size = args.NUM_IND)

    mutation_op = iaf.get_mutation_operator(m_rate = args.M_RATE,
                                            loci = 2,
                                            nrep = args.NUM_REP,
                                            burnin = args.burnin)

    output_op = iaf.get_output_operator(size = args.NUM_IND,
                                        m_rate = args.M_RATE,
                                        r_rate = args.R_RATE,
                                        s_rate = args.S_RATE,
                                        loci = 2,
                                        nrep = args.NUM_REP,
                                        ngen = args.NUM_GEN,
                                        burnin = args.burnin,
                                        output = args.OUTFILE)


    simulator = simu.Simulator(pops = pop)


    simulator.evolve(
        initOps = [init_info_op, init_genotype_op],
        preOps = mutation_op,
        matingScheme = mating_op,
        finalOps = output_op,
        gen = args.NUM_GEN + args.burnin)
