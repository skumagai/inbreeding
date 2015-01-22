"""
# selfingsim.utils

Small utility functions appearing in multiple modules.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

try:
    str = unicode
except NameError:
    pass

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

def getmode(mode):
    """
    Makes io.open, unicode, and csv work together in python2.

    Under this enviornment, a file must be opened in binary mode.
    AFAIU, io.open, when opened in text mode ("w"), expects unicode sequence,
    but csv module under python2 does not support unicode.
    A workaround is use python2's builtin open(), or use io.open in binary mode.
    I'm taking the latter approach.

    This function determing if the current interpreter is python2.  If so, It appends
    binary mode flag ("b") to the argument (e.g., "w" -> "wb", "a" -> "ab").

    """
    if sys.version_info.major < 3:
        return mode + "b"
    else:
        return mode
