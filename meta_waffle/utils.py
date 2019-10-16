"""
"""

import os
import errno
from datetime    import datetime
from time        import time


def mkdir(dnam):
    dnam = os.path.abspath(dnam)
    try:
        os.mkdir(dnam)
    except OSError as exc:
        if exc.errno != errno.EEXIST or not os.path.isdir(dnam):
            raise


def printime(msg, silent=True):
    if silent:
        return
    print (msg +
           (' ' * (79 - len(msg.replace('\n', '')))) +
           '[' +
           str(datetime.fromtimestamp(time()).strftime('%Y-%m-%d %H:%M:%S')) +
           ']')
