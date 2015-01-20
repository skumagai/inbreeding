from __future__ import print_function

# standard imports
import argparse
import json
import random
import sys

# within-package import
import data
import utils

def main():
    "Run this script as stand-alone."
    p = argparse.ArgumentParser()
    setup_command_line(sp)
    args = p.parse_args()
    args.func(args)

def setup_command_line(sp):
    p = sp.add_parser("sample")

    p.add_argument(
            "simfile",
            type = str,
            help = "file containing simulation results")
    p.add_argument(
            "generation",
            type = int,
            help = "sampling generation")
    p.add_argument(
            "samplesize",
            type = int,
            help = "sample size")
    p.add_argument(
            "reps",
            type = int,
            help = "number of replicates")
    p.set_defaults(func = sample)

    p = sp.add_parser("subsample")
    p.add_argument(
            "samplefile",
            type = str,
            help = "file containing sample (output of sample command)")
    p.add_argument(
            "subsamplefile",
            type = str,
            help = "file storing subsamples")
    p.add_argument(
            "samplesize",
            type = int,
            help = "sample size")
    p.add_argument(
            "sampleloci",
            type = int,
            nargs = "*",
            help = "indicies of loci to sample (0-based index)")
    p.set_defaults(func = subsample)

def sample(a):
    sims = data.createsample(a.simfile, a.generation)
    fbase = a.simfile.split(".")[:-1]

    nl = utils.getnewlinechar(a)

    nd = _ndigits(a.reps)
    ns = _ndigits(len(sims))

    template = fbase + ".simrep{{:0{}}}.{}.{{:0{}}}.json"

    for i, sim in enumerate(sims):
        for j in xrange(a.reps):
            s = sim.sample(a.samplesize)
            ofname = template.format(ns, a.samplesize, nd).format(i, j)
            with open(ofname, "w") as f:
                print(s.tojson(), sep = nl, end = nl, file = f)

def subsample(a):
    sample = data.createsample(a.samplefile)
    subsample = sample.sample(a.samplesize)

    nl = utils.getnewlinechar(a)

    with open(a.subsamplefile, "w") as f:
        print(subsample.tojson(), sep = nl, end = nl, file = f)

def _ndigits(n):
    digits = 1
    n = float(n)
    while n > 10.:
        n /= 10.
        digits += 1
    return digits

if __name__ == '__main__':
    main()
