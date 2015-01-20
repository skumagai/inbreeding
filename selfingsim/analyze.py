from __future__ import print_function

# standard imports
import argparse
from collections import Counter
import math

# within-package imports
import data
import utils

def main():
    sp = argparse.ArgumentParser().add_subparsers()
    setup_command_line(sp)
    a = sp.parse_args()
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
    sample = data.createsample(a.samplefile)
    src = sample.source

    if a.with_header:
        print("dataset\tsample\tselfing.gen")
    for i, t in zip(sample.ids, sample.tselfing):
        print("{}\tsample.{}\t{}".format(src, i, t))

if __name__ == '__main__':
    main()
