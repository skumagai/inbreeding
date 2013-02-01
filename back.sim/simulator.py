# -*- mode: python; coding: utf-8; -*-

# simulator.py - brief description

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

import driver as drv
import population as pop
import mutation as mut

class Simulator(object):
    """
    """

    def __init__(self, ploidy):
        """

        Arguments:
        - `ploidy`: ploidy of organisms
        """
        if 1 <= ploidy <= 2:
            self._driver = drv.createDriver(ploidy)
            self._population = pop.createPopulation(ploidy)
        else:
            raise Error


    def populationSize(self, size):
        """

        Arguments:
        - `size`: Population size at the time of sampling
        """
        self._population.size(size)

    def sampleSize(self, size, at=True, genes=1):
        """

        Arguments:
        - `size`: Number of sampled individuals
        - `at`: sampling deme
        - `genes': Number of sampled genes at a locus per sample
        """
        self._population.sampleSize = size


    def recombination(self, r):
        """

        Arguments:
        - `r`: probability of an offspring being recombinant
        """
        pass

    # def mutation(self, mu, at=True):
    #     """

    #     Arguments:
    #     - `mu`: per-site mutation rate per generation
    #     - `at`: applicable deme (True:= all demes)
    #     """
    #     pass

    def migration(self, m, src, dest=True):
        """

        Arguments:
        - `m`: migration rate
        - `src`: source
        - `dest`: destiation (True: to all other demes)
        """
        pass

    def registerEvent(self):
        """

        Arguments:
        - `self`:
        """
        pass

    def run(self, reps=1):
        """

        Arguments:
        - `reps`: number of replication
        """
        for n in range(reps):
            tree = self._driver.drive(pop)
            tree.print()
