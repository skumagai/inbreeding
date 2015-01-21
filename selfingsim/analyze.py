from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# standard imports
import argparse
from collections import Counter
import math

# within-package imports
import data
import utils

def main():
    p = argparse.ArgumentParser()
    sp = p.add_subparsers()
    setup_command_line(sp)
    a = p.parse_args()
    a.func(a)

def setup_command_line(sp):
    # setup shared command line argument
    pp = argparse.ArgumentParser(add_help = False)
    pp.add_argument(
            "samplefile",
            type = str,
            help = "file containing samples (output of (sub)sample command")
    pp.add_argument(
            "--with-header",
            action = "store_true",
            help = "set this flag to print headder line")

    # setup command line arguments for indiviudal subcommands
    p = sp.add_parser("inbcoeff", parents = [pp])
    p.set_defaults(func = inbcoeff)

    p = sp.add_parser("inbtime", parents = [pp])
    p.set_defaults(func = inbtime)

def inbcoeff(a):
    """
    Computes and prints inbreeding coefficients, Fis, and other related statistics
    for each locus as well as overall.

    Other statistics are observed and expected heterozygosities, bias-corrected
    inbreeding coefficients, and number of alleles.
    """
    sample = data.createsample(a.samplefile)
    coeffs = sample.inbreedingcoefficient()
    src = sample.source

    if a.with_header:
        print("dataset\tlocus\thetero.obs\thetero.exp\tFis\tFis.corrected\tnumber.of.alleles")
    for entry in coeffs:
        print("{}\t{}\t{}\t{}\t{}\t{}\t{}".
                format(
                    src,
                    entry["key"],
                    entry["hobs"],
                    entry["hexp"],
                    entry["fis"],
                    entry["fisc"],
                    entry["nalleles"]
                    )
                )

def inbtime(a):
    """
    Computes and prints number of generations untill the first outcrossing
    event along pedigree.

    This statistics are reported per-individual.
    If the most recent mating is outcrossing, this function returns 0.
    """
    sample = data.createsample(a.samplefile)
    src = sample.source

    if a.with_header:
        print("dataset\tsample\tselfing.gen")
    for i, t in zip(sample.ids, sample.tselfing):
        print("{}\tsample.{}\t{}".format(src, i, t))

if __name__ == '__main__':
    main()
