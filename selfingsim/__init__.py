"""
selfingsim
==========

Runs forward-in-time population genetic simulations for partially selfing organisms.
"""

import argparse

from . import analyze, convert, sample, simulate


def run():
    """
    Entry points of selfingsim
    """
    parser = argparse.ArgumentParser()
    subparsers = parser.add_subparsers()
    simulate.setup_command_line(subparsers)
    convert.setup_command_line(subparsers)
    sample.setup_command_line(subparsers)
    analyze.setup_command_line(subparsers)
    args = parser.parse_args()
    args.func(args)
