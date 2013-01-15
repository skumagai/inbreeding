# -*- mode: python; coding: utf-8; -*-

# driver.py - brief description

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

class Driver(object):
    """
    """

    def __init__(self, samples):
        """

        Arguments:
        - `samples`:
        """
        self._samples = samples

    def drive(self):
        """

        Arguments:
        - `self`:
        """
        samples = self._samples

        while samples.isCoalesced() is True:
            for ind in samples:
                if ind.isFromSelfing is True:
                    parent = ind.getParent()
                else:
                    if ind.numberOfInformativeChromosomes == 1:
                        parent = ind.getParent()
                    else:
                        (parent0, parent1) = ind.getParents()
            samples.incrementGeneration()


def main():
    """
    """
    pass

if __name__ == '__main__':
    main()
