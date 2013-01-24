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
import summarize

import mock

class TestFunctions(unittest.TestCase):

    def setUp(self):
        self.inds = [
            [([0],[1]), ([0],[1])],
            [([0],[1]), ([1],[0])],
            [([1],[0]), ([1],[0])]]


    def test_basic_info(self):
        self.assertEqual(summarize.basic_info(self.inds), (3, 2))


    def test_add_lists(self):
        lst0 = [0] * 5
        lst1 = [10] * 5
        summarize.add_lists(lst0, lst1)
        self.assertEqual(lst0, lst1)

        lst2 = [10] * 6
        self.assertRaises(IndexError, summarize.add_lists, lst0, lst2)


    def test_compute_f(self):
        self.assertEqual(summarize.compute_f(self.inds), [2./3., 2./3.])


    def test_compute_g(self):
        self.assertEqual(summarize.compute_g(self.inds), [1./3., 1./3.])


    def test_compute_P(self):
        exp = {
            (0,0): 1./3.,
            (1,1): 2./3.,
            (0,1): 0,
            (1,0): 0.
            }
        self.assertEqual(summarize.compute_P(self.inds), exp)


    def test_compute_W(self):
        exp = {
            (0,0): 8./12.,
            (1,1): 4./12.,
            (0,1): 0.,
            (1,0): 0.
            }
        self.assertEqual(summarize.compute_W(self.inds), exp)


    def test_check_identity(self):
        self.assertEqual(summarize.check_identity(2, self.inds[0][1], self.inds[1][1]), [0, 0])
