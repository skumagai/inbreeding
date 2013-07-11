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


import simuPOP as simu
import csv

def get_population(size, loci, infoFields='self_gen'):
    """Construct a population object."""
    return simu.Population(size = size,
                           ploidy = 2,
                           loci = loci,
                           infoFields = infoFields)


def get_init_info(field='self_gen'):
    """Zero initialize info field `field`."""
    return simu.InitInfo(0, infoFields=field)


def get_init_genotype():
    """Zero initialize genotype of all organisms."""
    return simu.InitGenotype(prop=[1])


def pickTwoParents(pop):
    """
    Pick two distinct individuals to be parents for outcrossing.
    """
    parents = list(pop.individuals())
    n_par = len(parents)
    rng = simu.getRNG()
    while True:
        par0 = rng.randInt(n_par)
        par1 = rng.randInt(n_par)
        while par0 == par1:
            par1 = rng.randInt(n_par)
        yield (par0, par1)


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
        Increment the value of offspring's infofield `field` by one.

        This method does not use only `dad` not `mom`.
        """
        off.setInfo(dad.info(self.field) + 1, self.field)
        return True


class MyOutcrossingTagger(simu.PyOperator):
    """
    Update information field to reflex outcrossing.

    When outcrossing occurred, this operator reset the value of `field` to indicate
    that an offspring was from outcrossing rather than selfing.
    """

    def __init__(self, field='self_gen'):
        self.field = field

        super(MyOutcrossingTagger, self).__init__(func = self.record)

    def record(self, pop, off, dad, mom):
        """
        Reset the value of offspring's infofield `field` to zero.

        This method does not use `dad` and `mom`.
        """
        off.setInfo(0, self.field)
        return True
