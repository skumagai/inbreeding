"""
selfingsim.convert
==================

Conversion between various file formats.
"""


# standard imports
import argparse
import io

# within-package imports
from . import data, utils


def run():
    """
    Runs file format conversion as a stand-alone script.
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
    # setup shared command line arguments
    sharedparser1 = argparse.ArgumentParser(add_help=False)
    sharedparser1.add_argument(
        "samplefile",
        type=str,
        help="file containing samples (output of (sub)sample command)",
    )

    sharedparser2 = argparse.ArgumentParser(add_help=False)
    sharedparser2.add_argument(
        "-w",
        action="store_true",
        default=False,
        help="set this flag to force Windows newline character (CRLF).",
    )

    # setup command line arguments for individual subcommands
    subparser = subparsers.add_parser("phase", parents=[sharedparser1, sharedparser2])
    subparser.set_defaults(func=phase)

    subparser = subparsers.add_parser("phase2rmes", parents=[sharedparser2])
    subparser.add_argument(
        "phasefile", type=str, help="file containing samples in phase format."
    )
    subparser.set_defaults(func=phase2rmes)

    subparser = subparsers.add_parser("nexus", parents=[sharedparser1, sharedparser2])
    subparser.add_argument(
        "localsizes", type=int, nargs="*", help="local sample sizes (zero or more)"
    )
    subparser.set_defaults(func=nexus)

    subparser = subparsers.add_parser("rmes", parents=[sharedparser1, sharedparser2])
    subparser.set_defaults(func=rmes)

    subparser = subparsers.add_parser("rmescombine", parents=[sharedparser2])
    subparser.add_argument("combinedfile", type=str, help="file of output")
    subparser.add_argument(
        "rmesfiles", type=str, nargs="+", help="one or more sample files in RMES format"
    )
    subparser.set_defaults(func=rmescombine)


# The following functions are mainly thin-wrapper of formatting methods in
# data.{Basic,Full}Sample objects.


def phase(config):
    """
    Creates a PHASE-formatted file.

    The result can be used as an input to BALI-PHY.
    """
    ofname = ".".join(config.samplefile.split(".")[:-1] + ["phase"])
    sample = data.createsample(config.samplefile)
    newline = utils.getnewlinechar(config)
    with io.open(ofname, "w") as fhandle:
        print(*sample.tophase(), sep=newline, end=newline, file=fhandle)


def nexus(config):
    """
    Creates a NEXUS-formatted file.

    The result can be used as an input to GDA.
    """
    ofname = ".".join(config.samplefile.split(".")[:-1] + ["nex"])
    sample = data.createsample(config.samplefile)
    newline = utils.getnewlinechar(config)
    with io.open(ofname, "w") as fhandle:
        print(*(sample.tonexus()), sep=newline, end=newline, file=fhandle)


def rmes(config):
    """
    Creates a RMES-formatted file.

    The result can be used as an input to RMES.
    """
    ofname = ".".join(config.samplefile.split(".")[:-1] + ["rmes"])
    sample = data.createsample(config.samplefile)
    newline = utils.getnewlinechar(config)
    with io.open(ofname, "w") as fhandle:
        print(*(sample.tormes()), sep=newline, end=newline, file=fhandle)


def rmescombine(config):
    """
    Combines multiple RMES-formatted files.

    The result can still be used as an input to RMES.
    """
    samples = data.tormes([data.createsample(fname) for fname in config.rmesfiles])
    newline = utils.getnewlinechar(config)
    with io.open(config.combinedfile, "w") as fhandle:
        print(*samples, sep=newline, end=newline, file=fhandle)


def phase2rmes(config):
    """
    Converts a PHASE-formatted file to a RMES-formatted file.
    """
    ofname = ".".join(config.samplefile.split(".")[:-1] + ["rmes"])
    sample = data.createsample(config.samplefile)
    newline = utils.getnewlinechar(config)
    with io.open(ofname, "w") as fhandle:
        print(*(sample.tormes()), sep=newline, end=newline, file=fhandle)


if __name__ == "__main__":
    run()
