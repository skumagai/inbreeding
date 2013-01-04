# -*- mode: python; coding: utf-8; -*-

# utility.py - brief description

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

from __future__ import print_function

import json
import os.path


def import_right_module(args):
    mode = 'infinite-sites'
    try:
        global sim
        import simuOpt
        info = get_info(args.DIR[0])
        if u'mode' in info and info[u'mode'] == u'infinite-alleles':
            mode = 'infinite-alleles'
            simuOpt.setOptions(alleleType = 'long')
        else:
            simuOpt.setOptions(alleleType = 'binary')
        import simuPOP as sim
    except ImportError as e:
        print('[ERROR] {}'.format(e), file=sys.stderr)
        sys.exit(1)
    return mode


def get_info(dir):
    with open(os.path.join(dir, 'conf.json'), 'r') as rf:
        info = json.load(rf)
    return info


def chunks(l, n):
    '''Divide a list into equal-length sublists'''
    for i in xrange(0, len(l), n):
        yield l[i:i+n]
