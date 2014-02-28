# -*- mode: python; coding: utf-8; -*-

# main.py - An entry point to partial selfing simulations.  The
# simulations are conducted with population genetics forward simulator
# simuPOP, and four scenarios are explored.

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

import sys
import argparse


# Selectively import simuPOP with an appropriate alleleType.  Because
# python caches imported modules, subsequent import of simuPOP in
# model specific submodules is not re-imported.  Therefore, those
# submodules will automatically use the right version of simuPOP.
def exec_infinite_sites(args):
    import partial_selfing.infinite_sites as model
    model.run(args)


# See the comment in front of exec_two_loci
def exec_infinite_alleles(args):
    import partial_selfing.infinite_alleles as model
    model.run(args)


def summarize(args):
    import partial_selfing.summary as model
    model.run(args)


if __name__ == '__main__':

    # Handling model specification through command line arguments.
    #
    # No argument has a default value, so all values have to be
    # explicitly supplied.  This is by design.
    #
    # Depending on which scenario, 2 loci with partial linkage or many
    # loci with free linkage, different functions are invoked.  Those
    # functions import an approximate version of simuPOP module.

    parser_common = argparse.ArgumentParser(add_help=False)
    parser_common.add_argument('OUTFILE',
                               type=str,
                               help='path to output file')
    parser_common.add_argument('NUM_IND',
                               type=int,
                               help='number of individuals')
    parser_common.add_argument('NUM_GEN',
                               type=int,
                               help='number of generations')
    parser_common.add_argument('NUM_REP',
                               type=int,
                               help='number of replicates')
    parser_common.add_argument('NUM_LOCI',
                               type=int,
                               help='number of loci')
    parser_common.add_argument('-t',
                               dest='M_RATE',
                               required=True,
                               nargs='+',
                               type=float,
                               help='mutation rate')
    parser_common.add_argument('-s',
                               dest='S_RATE',
                               required=True,
                               type=float,
                               help='selfing rate')
    parser_common.add_argument('-r',
                               dest='R_RATE',
                               required=True,
                               type=float,
                               help='recombination rate')
    parser_common.add_argument('--burnin',
                               type=int,
                               default=0,
                               help='length of burnin phase')
    parser_common.add_argument('--debug',
                               type=int,
                               default=0,
                               help='debug output')
    parser_common.add_argument('--output_per',
                               type=int,
                               default=0,
                               help='frequency of outputting population state')

    parser = argparse.ArgumentParser(description='Foward simulate partial selfing')

    subparsers = parser.add_subparsers()

    parser_alleles = subparsers.add_parser('infinite_alleles',
                                           help='infinite-alleles model',
                                           parents=[parser_common])
    parser_alleles.set_defaults(func=exec_infinite_alleles)

    parser_alleles.add_argument('--distinct_init',
                                type=int,
                                default=0,
                                help='initial population is not monomorphic')

    parser_sites = subparsers.add_parser('infinite_sites',
                                         help='infinite-sites model',
                                         parents=[parser_common])
    parser_sites.add_argument('--allele_length',
                              type=int,
                              default=2**6, # 64
                              help='number of maximum polymorphic sites at one time')
    parser_sites.set_defaults(func=exec_infinite_sites)


    summary = subparsers.add_parser('summary',
                                    help='summarize simulation results')

    summary.add_argument('PATTERN',
                         type=str,
                         help='glob pattern of simulation result files')
    summary.set_defaults(func=summarize)


    args = parser.parse_args()
    if len(args.M_RATE) == 1:
        args.M_RATE = [args.M_RATE[0] for i in range(args.NUM_LOCI)]
    elif len(args.M_RATE) != args.NUM_LOCI:
        sys.stderr.write("Number of mutation rates must be equal to 1 or the number of loci.\n")
        parser.print_help()
        sys.exit(2)

    # Enter into an appropriate entry point.
    args.func(args)
