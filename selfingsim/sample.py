"""
selfingsim.sample
=================

Sample subsets of organisms from simulation results or other samples.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# standard imports
import argparse
import io

# within-package import
from . import data
from . import utils

def run():
    """
    Runs this script as stand-alone.
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    setup_command_line(subparsers)
    args = parser.parse_args()
    args.func(args)

def setup_command_line(subparsers):
    """
    Sets up command line arguments.
    """
    parser = subparsers.add_parser("sample")
    parser.add_argument(
        "simfile",
        type=str,
        help="file containing simulation results")
    parser.add_argument(
        "generation",
        type=int,
        help="sampling generation")
    parser.add_argument(
        "samplesize",
        type=int,
        help="sample size")
    parser.add_argument(
        "reps",
        type=int,
        help="number of replicates")
    parser.set_defaults(func=sample)

    parser = subparsers.add_parser("subsample")
    parser.add_argument(
        "samplefile",
        type=str,
        help="file containing sample (output of sample command)")
    parser.add_argument(
        "subsamplefile",
        type=str,
        help="file storing subsamples")
    parser.add_argument(
        "samplesize",
        type=int,
        help="sample size")
    parser.add_argument(
        "sampleloci",
        type=int,
        nargs="*",
        help="indicies of loci to sample (0-based index)")
    parser.set_defaults(func=subsample)

def sample(config):
    """
    Samples a subset of individuals from the current sample from simulation
    results in TSV format.
    """
    # Assume that one simulation result does not contain results of more than one replicates.
    sim = data.createsample(config.simfile, config.generation)[0]
    fbase = ".".join(config.simfile.split(".")[:-1])

    dreps = _ndigits(config.reps)

    # the following template is used for having the right amount of padding
    # in the output file name.
    template = fbase + ".size_{}.sample_rep_{{:0{}}}.json"

    for j in xrange(config.reps):
        newsamp = sim.sample(config.samplesize)
        ofname = template.format(config.samplesize, dreps).format(j)
        with io.open(ofname, "w") as fhandle:
            print(newsamp.tojson(), file=fhandle)

def subsample(config):
    """
    Gets a subsample from already a sample in a JSON-formatted file.
    """
    samp = data.createsample(config.samplefile)
    subs = samp.sample(config.samplesize, config.sampleloci)

    with io.open(config.subsamplefile, "w") as fhandle:
        print(subs.tojson(), file=fhandle)

def _ndigits(number):
    """
    Count the number of digits required to represent in decimal.
    """
    digits = 1
    while number > 10.:
        number /= 10.
        digits += 1
    return digits

if __name__ == '__main__':
    run()
