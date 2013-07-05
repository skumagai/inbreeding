# -*- mode: python; coding: utf-8; -*-

# pself.py - An entry point to partial selfing simulations.  The
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

import argparse


# Selectively import simuPOP with an appropriate alleleType.  Because
# python caches imported modules, subsequent import of simuPOP in
# model specific submodules is not re-imported.  Therefore, those
# submodules will automatically use the right version of simuPOP.
def exec_two_loci(args):
    if args.M_TYPE == 'infinite_loci':
        import pself.two_loci.infinite_sites as model
    else:                       # infinite_alleles
        import pself.two_loci.infinite_alleles as model
    model.run(args)


# See the comment in front of exec_two_loci
def exec_many_loci(args):
    if args.M_TYPE == 'infinite_loci':
        import pself.many_loci.infinite_sites as model
    else:                       # infinite_alleles
        import pself.many_loci.infinite_alleles as model
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
    parser_common.add_argument('NUM_IND',
                               type=int,
                               help='number of individuals')
    parser_common.add_argument('NUM_GEN',
                               type=int,
                               help='number of generations')
    parser_common.add_argument('NUM_REP',
                               type=int,
                               help='number of replicates')
    parser_common.add_argument('M_TYPE',
                               choices=['infinite_loci', 'infinite_alleles'],
                               help='mutational model')
    parser_common.add_argument('M_RATE',
                               type=float,
                               help='mutation rate')
    parser_common.add_argument('--burnin',
                               type=int,
                               default=0,
                               help='length of burnin phase')

    parser = argparse.ArgumentParser(description='Foward simulate partial selfing')

    subparsers = parser.add_subparsers()

    parser_two = subparsers.add_parser('two_loci',
                                       help='two partially linked loci',
                                       parents=[parser_common])
    parser_two.add_argument('R_RATE', type=float, help='recombination rate')
    parser_two.set_defaults(func=exec_two_loci)

    parser_many = subparsers.add_parser('many_loci',
                                        help='many unliked loci',
                                        parents=[parser_common])
    parser_many.add_argument('NUM_LOCI',
                             type=int,
                             help='number of loci')
    parser_many.add_argument('--allele_length',
                             type=int,
                             defualt=2**6, # 64
                             help='number of maximum polymorphic sites at one time')
    parser_many.set_defaults(func=exec_many_loci)

    args = parser.parse_args()

    # Enter into an appropriate entry point.
    args.func(args)
