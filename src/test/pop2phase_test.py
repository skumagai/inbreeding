# -*- mode: python; coding: utf-8; -*-

# pop2phase_test.py - brief description

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

'''unittest for pop2phase.py.  Some top-level or nearly-so functions are not tested'''


from __future__ import print_function

import unittest
import random
import argparse
import itertools
import contextlib
import simuPOP as sim
import csv
import pop2phase

import mock

class TestFunctions(unittest.TestCase):

    def setUp(self):
        self.args = argparse.Namespace()
        self.args.DIR = 'stub'


    def test_get_mode(self):
        with mock.patch('pop2phase.get_info') as info_mock:
            info_mock.return_value = {u'mode': u'infinite-sites'}
            self.assertEqual(pop2phase.get_mode(self.args), 'infinite-sites')

        with mock.patch('pop2phase.get_info') as info_mock:
            info_mock.return_value = {u'mode': u'infinite-alleles'}
            self.assertEqual(pop2phase.get_mode(self.args), 'infinite-alleles')

        with mock.patch('pop2phase.get_info') as info_mock:
            info_mock.return_value = {u'mode': None}
            self.assertEqual(pop2phase.get_mode(self.args), 'infinite-sites')


    def test_encode(self):
        codes = {'idx': 0}
        self.assertEqual(0, pop2phase.encode([0,0,0], codes))
        self.assertEqual(1, pop2phase.encode([0,1,0], codes))
        self.assertEqual(2, pop2phase.encode([1,1,0], codes))
        self.assertEqual(1, pop2phase.encode([0,1,0], codes))
        self.assertEqual(0, pop2phase.encode([0,0,0], codes))
        self.assertRaises(ValueError, pop2phase.encode, [0,0], codes)


    def test_convert_genotype(self):
        loci_dict = {i: {'idx': 0} for i in xrange(3)}
        g0 = pop2phase.convert_genotype([[0,1], [2,3], [4,5]], loci_dict)
        g1 = pop2phase.convert_genotype([[0,1], [2,3], [4,5]], loci_dict)
        self.assertEqual(list(itertools.chain.from_iterable(zip(g0, g1))),
                         [0 for i in xrange(6)])
        g0 = pop2phase.convert_genotype([[0,1], [2,3], [4,5]], loci_dict)
        g1 = pop2phase.convert_genotype([[1,1], [2,3], [5,5]], loci_dict)
        self.assertEqual(list(itertools.chain.from_iterable(zip(g0, g1))),
                         [0, 1, 0, 0, 0, 1])
        g0 = pop2phase.convert_genotype([[1,0], [1,3], [6,5]], loci_dict)
        g1 = pop2phase.convert_genotype([[0,1], [2,3], [4,5]], loci_dict)
        self.assertEqual(list(itertools.chain.from_iterable(zip(g0, g1))),
                         [2, 0, 1, 0, 2, 0])



    def test_sample_loci(self):

        with  mock.patch('__builtin__.open', create=True) as mock_open, \
          mock.patch('csv.reader', create=True) as mock_csv_reader, \
          mock.patch('random.sample') as sample_mock:
            mock_open.return_value = mock.MagicMock(spec=file)
            mock_csv_reader.return_value = mock.MagicMock()
            mock_csv_reader.return_value.next.return_value = ['5']
            mock_csv_reader.return_value.__iter__.return_value = \
              iter([list('0101010101'),
                    list('1010101010'),
                    list('0121012101'),
                    list('2101210121'),
                    list('0101010101')])
            sample_mock.return_value = [1,3]

            self.assertEqual(pop2phase.sample_loci('stub', 2, 2),
                             [('1', '0', '1', '0'),
                              ('0', '1', '0', '1')])


    def test_randomize(self):
        pass
