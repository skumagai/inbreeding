"""
# selfingsim.utils

Small utility functions appearing in multiple modules.
"""


import os
import sys


def getnewlinechar(config):
    """
    Returns a new line character.
    """
    if config.w:
        return "\r\n"
    else:
        return str(os.linesep)
