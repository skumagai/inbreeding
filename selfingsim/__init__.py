"""
selfingsim
==========

Runs forward-in-time population genetic simulations for partially selfing organisms.
"""
from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import argparse

from . import simulate
from . import convert
from . import sample
from . import analyze

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

