"""
selfingsim.analyze
##################

Analyze simulation data.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# standard imports
import argparse

# within-package imports
from . import data

def run():
    """
    Runs analysis as a stand-alone script.
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    setup_command_line(subparsers)
    args = parser.parse_args()
    args.func(args)

def setup_command_line(subparsers):
    """
    Sets up command line interface.
    """
    # setup shared command line argument
    sharedparser = argparse.ArgumentParser(add_help=False)
    sharedparser.add_argument(
        "samplefile",
        type=str,
        help="file containing samples (output of (sub)sample command")
    sharedparser.add_argument(
        "--with-header",
        action="store_true",
        help="set this flag to print headder line")

    # setup command line arguments for indiviudal subcommands
    subparser = subparsers.add_parser("inbcoeff", parents=[sharedparser])
    subparser.set_defaults(func=inbcoeff)

    subparser = subparsers.add_parser("inbtime", parents=[sharedparser])
    subparser.set_defaults(func=inbtime)

def inbcoeff(config):
    """
    Computes and prints inbreeding coefficients, Fis, and other related statistics
    for each locus as well as overall.

    Other statistics are observed and expected heterozygosities, bias-corrected
    inbreeding coefficients, and number of alleles.
    """
    sample = data.createsample(config.samplefile)
    coeffs = sample.inbreedingcoefficient()
    src = sample.source

    if config.with_header:
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
                  entry["nalleles"]))

def inbtime(config):
    """
    Computes and prints number of generations untill the first outcrossing
    event along pedigree.

    This statistics are reported per-individual.
    If the most recent mating is outcrossing, this function returns 0.
    """
    sample = data.createsample(config.samplefile)
    src = sample.source

    if config.with_header:
        print("dataset\tsample\tselfing.gen")
    for i, tselfing in zip(sample.ids, sample.tselfing):
        print("{}\tsample.{}\t{}".format(src, i, tselfing))

if __name__ == '__main__':
    run()
