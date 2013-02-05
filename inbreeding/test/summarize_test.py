# -*- mode: python; coding: utf-8; -*-

# summarize_test.py - brief description

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

'''unittest for summarize.py.  Some top-level or nearly-so functions are not tested'''

import unittest
from inbreeding import summarize

import mock

class TestFunctions(unittest.TestCase):

    def setUp(self):
        self.inds = {(((0,),(1,)), ((0,),(1,))): 2,
                     (((0,),(1,)), ((1,),(0,))): 3,
                     (((1,),(0,)), ((1,),(0,))): 4}
        self.num_loci = 2
        self.num_individuals = float(9)


    def test_add_lists(self):
        lst0 = [0] * 5
        lst1 = [10] * 5
        self.assertNotEqual(lst0, lst1)
        summarize.add_lists(lst0, lst1)
        self.assertEqual(lst0, lst1)

        lst2 = [10] * 6
        self.assertRaises(IndexError, summarize.add_lists, lst0, lst2)


    def test_compute_f(self):
        self.assertEqual(summarize.compute_f(self.num_loci, self.inds), [6, 6])


    def test_compute_g(self):
        self.assertEqual(summarize.compute_g(self.num_loci, self.inds), [70, 70])


    def test_compute_P(self):
        exp = {
            (0,0): 6,
            (1,1): 3,
            (0,1): 0,
            (1,0): 0
            }
        self.assertEqual(summarize.compute_P(self.num_loci, self.inds), exp)


    def test_compute_W(self):
        exp = {
            (0,0): 70,
            (1,1): 74,
            (0,1): 0,
            (1,0): 0
            }
        self.assertEqual(summarize.compute_W(self.num_loci, self.inds), exp)


    def test_check_identity(self):
        keys = sorted(self.inds.keys())
        self.assertEqual(summarize.check_identity(2, keys[0][0], keys[0][1]), [0, 0])
        self.assertEqual(summarize.check_identity(2, keys[0][1], keys[1][1]), [1, 1])
        self.assertEqual(summarize.check_identity(2, ((1,),(0,)), ((0,), (0,))), [1, 0])


    def test_relationship_f_P(self):
        f = summarize.compute_f(self.num_loci, self.inds)
        f = [v / self.num_individuals for v in f]
        P = summarize.compute_P(self.num_loci, self.inds)
        total = float(sum(val for val in P.values()))
        P = {i: v / total for i, v in P.items()}
        self.assertAlmostEqual(f[0], P[(0,0)] + P[(0,1)])
        self.assertAlmostEqual(1 - f[0], P[(1,0)] + P[(1,1)])

        self.assertAlmostEqual(f[1], P[(0,0)] + P[(1,0)])
        self.assertAlmostEqual(1 - f[1], P[(0,1)] + P[(1,1)])

    def test_relationship_g_W(self):
        g = summarize.compute_g(self.num_loci, self.inds)
        total = self.num_individuals * (self.num_individuals - 1) * 2
        g = [i / total for i in g]

        W = summarize.compute_W(self.num_loci, self.inds)
        total = float(sum(val for val in W.values()))
        W = {i: v / total for i, v in W.items()}
        self.assertAlmostEqual(g[0], W[(0,0)] + W[(0,1)])
        self.assertAlmostEqual(1 - g[0], W[(1,0)] + W[(1,1)])

        self.assertAlmostEqual(g[1], W[(0,0)] + W[(1,0)])
        self.assertAlmostEqual(1 - g[1], W[(0,1)] + W[(1,1)])
