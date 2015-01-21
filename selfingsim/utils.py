from __future__ import absolute_import
from __future__ import division
from __future__ import print_function
from __future__ import unicode_literals

def getnewlinechar(a):
    if a.w:
        return "\r\n"
    else:
        return os.linesep
