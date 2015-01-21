from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

# standard imports
import argparse
import io
import itertools
import os
import os.path

# within-package imports
import data
import utils

def main():
    s = argparse.ArgumentParser()
    sp = s.add_subparsers()
    setup_commnad_line(sp)
    a = p.parse_args()
    a.func(a)

def setup_command_line(sp):
    # setup shared command line arguments
    pp1 = argparse.ArgumentParser(add_help = False)
    pp1.add_argument(
            "samplefile",
            type = str,
            help = "file containing samples (output of (sub)sample command)")

    pp2 = argparse.ArgumentParser(add_help = False)
    pp2.add_argument(
            "-w",
            action = "store_true",
            default = False,
            help = "set this flag to force Windows newline character (CRLF).")

    # setup command line arguments for individual subcommands
    p = sp.add_parser("phase", parents = [pp1, pp2])
    p.set_defaults(func = phase)

    p = sp.add_parser("phase2rmes", parents = [pp2])
    p.add_argument(
            "phasefile",
            type = str,
            help = "file containing samples in phase format.")
    p.set_defaults(func = phase2rmes)

    p = sp.add_parser("nexus", parents = [pp1, pp2])
    p.add_argument(
            "localsizes",
            type = int,
            nargs = "*",
            help = "local sample sizes (zero or more)"
            )
    p.set_defaults(func = nexus)

    p = sp.add_parser("rmes", parents = [pp1, pp2])
    p.set_defaults(func = rmes)

    p = sp.add_parser("rmescombine", parents = [pp2])
    p.add_argument(
            "combinedfile",
            type =str,
            help = "file of output")
    p.add_argument(
            "rmesfiles",
            type = str,
            nargs = "+",
            help = "one or more sample files in RMES format")
    p.set_defaults(func = rmescombine)

# The following functions are mainly thin-wrapper of formatting methods in
# data.{Basic,Full}Sample objects.

def phase(a):
    """
    Creates a PHASE-formatted file.

    The result can be used as an input to BALI-PHY.
    """
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["phase"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with io.open(ofname, "w") as f:
        print(sample.tophase(), sep = nl, end = nl, file = f)

def nexus(a):
    """
    Creates a NEXUS-formatted file.

    The result can be used as an input to GDA.
    """
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["nex"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with io.open(ofname, "w") as f:
        print(*(sample.tonexus()), sep = nl, end = nl, file = f)

def rmes(a):
    """
    Creates a RMES-formatted file.

    The result can be used as an input to RMES.
    """
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["rmes"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with io.open(ofname, "w") as f:
        print(*(sample.tormes()), sep = nl, end = nl, file = f)

def rmescombine(a):
    """
    Combines multiple RMES-formatted files.

    The result can still be used as an input to RMES.
    """
    samples = [data.createsample(fname) for fname in a.rmesfiles]
    nl = utils.getnewlinechar(a)
    with io.open(a.combinedfile, "w") as f:
        print(
                *([l for sample in samples for line in sample.tormes()]),
                sep = nl,
                end = nl,
                file = f)

def phase2rmes(a):
    """
    Converts a PHASE-formatted file to a RMES-formatted file.
    """
    # TODO: Handle missing value (represented as "-9").
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["rmes"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with io.open(ofname, "w") as f:
        print(*(sample.tormes()), sep = nl, end = nl, file = f)

if __name__ == '__main__':
    main()
