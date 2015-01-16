# -*- mode: python; coding: utf-8; -*-

# post_process.py - collection of post process routines on simulated data

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

import csv, random
import os.path as path


def get_samples(N, size):
    return random.sample(xrange(N), size)


def phase(config, gen, i, samples):
    nsam = len(samples)
    nloci = len(samples[0][1])

    outfile = config.outfile
    outfile = path.splitext(outfile)[0] + 'gen{}.rep{}.phase.txt'.format(gen, i)

    with open(outfile, 'w') as ofp:
        ofp.write('{}\n{}\n'.format(nsam, nloci))
        ofp.write('{}\n'.format('M'*nloci))

        for s in samples:
            ofp.write('sample_{}\t{}\n'.format(s[0],
                                               '\t'.join('{}\t{}'.format(*i) for i in s[1])))



def rmes(config, gen, i, samples):
    nsam = len(samples)
    nloci = len(samples[0][1])

    outfile = config.outfile
    outfile = path.splitext(outfile)[0] + 'gen{}.rep{}.rmes.txt'.format(gen, i)

    with open(outfile, 'w') as ofp:
        ofp.write('{}\r\n{}\r\n{}\r\n{}\r\n'.format(outfile, nsam, nloci))
        for s in samples:
            ofp.write('{}\r\n'.format(' '.join(str(int(i[0] != i[1])) for i in s[1])))



def dispatch(config, command, sample):
    if command == 'phase':
        phase(config, sample)
    elif command == 'rmes':
        rmes(config, sample)
    else:
        print('undefined command: {}'.format(command))


def get_pops(config, gens):
    pops = dict()
    with open(config.outfile, 'r') as fp:
        reader = csv.reader(fp)

        hs = next(reader)
        gidx = hs.index('generation')
        l0idx = hs.index('locus 0')

        for r in next(reader):
            nr = next(reader)
            gen = r[gidx]
            if gen in gens:
                ind = list([(i, j) for i, j in zip(r[l0idx:], nr[l0idx:])])
                try:
                    pops[gen].append(ind)
                except:
                    pops[gen] = [ind]
        return pops


def run(config):
    if config.sample_at == 'last':
        gens = [str(config.gens)]
    else:
        gens = list(str(i) for i in range(0, config.gens + 1, config.output_per))
    pops = get_pops(config, gens)

    for pop, inds in pops.items():
        for i in range(config.reps):
            sidx = get_samples(config.N, config.sample_size)
            samples = [(i, inds[i]) for i in sidx]

            try:
                for command in config.command:
                    dispatch(config, command, samples)
            except:
                try:
                    dispatch(config, config.command, samples)
                except:
                    print('wrong argument for "command"')
