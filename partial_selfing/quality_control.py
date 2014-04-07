# -*- mode: python; coding: utf-8; -*-

# post_process.py - functions for post process simulation results

# Copyright (C) 2014 Seiji Kumagai

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

import ast, csv

GENROW = 11
G0 = 15

def to_prop(data):
    total = float(sum(i for i in data[0].values()))
    return {k: {kk: v / total for kk, v in data[k].items()}
            for k in data.keys()}

def spectra(args):
    path = args.path
    with open(path, 'r') as f, open(path.split('.')[0] +'.spec', 'w') as of:
        r = csv.reader(f)
        n = len(next(r)[G0:])
        gen = None

        for row in r:
            if gen != row[GENROW]:
                gen = row[GENROW]
                try:
                    of.write(str(to_prop(spec)) + '\n')
                    spec.clear()
                except NameError:
                    pass
                spec = {i: {} for i in range(n)}
            for i, col in enumerate(row[G0:]):
                try:
                    spec[i][int(col)] += 1
                except KeyError:
                    spec[i][int(col)] = 1
        of.write(str(to_prop(spec)) + '\n')


def split(args):
    path = args.path
    with open(path, 'r') as f:
        r = csv.reader(f)

        gen = None

        header = next(r)

        for row in r:
            if row[GENROW] != gen:
                gen = row[GENROW]
                try:
                    of.close()
                except NameError:
                    pass
                of = open(path.split('.')[0] + '.' + gen + '.csv', 'w')
                w = csv.writer(of)
                w.writerow(header)
            w.writerow(row)
        of.close()



def count_nalleles(args):
    path = args.path
    with open(path, 'r') as f, open(path.split('.')[0] + '.nalleles.csv', 'w') as of:
        w = csv.writer(of)
        for row in f:
            data = []
            for v in ast.literal_eval(row).values():
                data.append(len(v))
            w.writerow(data)

def count_zeros(args):
    path = args.path
    with open(path, 'r') as f, open(path.split('.')[0] + '.zeros.csv', 'w') as of:
        w = csv.writer(of)
        for row in f:
            data = []
            for v in ast.literal_eval(row).values():
                try:
                    data.append(v[0])
                except KeyError:
                    data.append(0.0)
            w.writerow(data)

def count_inits(args):
    path = args.path
    max_init = 2 * args.pop_size
    with open(path, 'r') as f, open(path.split('.')[0] + '.inits.csv', 'w') as of:
        w = csv.writer(of)
        for row in f:
            data = []
            for v in ast.literal_eval(row).values():
                frac = 0.0
                for i in range(max_init):
                    try:
                        frac += v[i]
                    except KeyError:
                        pass
                data.append(frac)
            w.writerow(data)

if __name__ == '__main__':

    import argparse

    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    spe = sp.add_parser('spectra')
    spl = sp.add_parser('split')
    cn = sp.add_parser('alleles')
    cz = sp.add_parser('zeros')
    ci = sp.add_parser('inits')

    spe.add_argument('path', type=str)
    spl.add_argument('path', type=str)
    cn.add_argument('path', type=str)
    cz.add_argument('path', type=str)
    ci.add_argument('path', type=str)
    ci.add_argument('pop_size', type=int)

    spe.set_defaults(func=spectra)
    spl.set_defaults(func=split)
    cn.set_defaults(func=count_nalleles)
    cz.set_defaults(func=count_zeros)
    ci.set_defaults(func=count_inits)

    args = p.parse_args()
    args.func(args)
