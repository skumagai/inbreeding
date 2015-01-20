from __future__ import print_function
# standard imports
import argparse
import itertools
import os
import os.path

# within-package imports
import data
import utils

def main():
    sp = argparse.ArgumentParser().add_subparsers()
    setup_commnad_line(sp)
    a = sp.parse_args()
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

# The following functions are mainly thing-wrapper of formatting methods in
# data.{Basic,Full}Sample objects.

def phase(a):
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["phase"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with open(ofname, "w") as f:
        print(sample.tophase(), sep = nl, end = nl, file = f)

def nexus(a):
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["nex"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with open(ofname, "w") as f:
        print(*(sample.tonexus()), sep = nl, end = nl, file = f)

def rmes(a):
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["rmes"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with open(ofname, "w") as f:
        print(*(sample.tormes()), sep = nl, end = nl, file = f)

def rmescombine(a):
    samples = [data.createsample(fname) for fname in a.rmesfiles]
    nl = utils.getnewlinechar(a)
    with open(a.combinedfile, "w") as f:
        print(
                *([l for sample in samples for line in sample.tormes()]),
                sep = nl,
                end = nl,
                file = f)

def phase2rmes(a):
    # TODO: Handle missing value (represented as "-9").
    ofname = ".".join(a.samplefile.split(".")[:-1] + ["rmes"])
    sample = data.createsample(a.samplefile)
    nl = utils.getnewlinechar(a)
    with open(ofname, "w") as f:
        print(*(sample.tormes()), sep = nl, end = nl, file = f)

