"""
# selfingsim.utils

Small utility functions appearing in multiple modules.
"""

from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

import os

def getnewlinechar(config):
    """
    Returns a new line character.
    """
    if config.w:
        return "\r\n"
    else:
        return unicode(os.linesep)
