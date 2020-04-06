"""Misc
"""

# Already in python_startup.py
import sys
import os

# Others
import time
# import re
# import builtins
import dataclasses as dcs
from collections import namedtuple
# from typing import Optional, Any


import time
class Timer():
    """Timer.

    Example::

      with Timer('<description>'):
        do_stuff()
    """
    def __init__(self, name=None):
        self.name = name

    def __enter__(self):
        self.tstart = time.time()

    def __exit__(self, type, value, traceback):
        #pass # Turn off timer messages
        if self.name:
            print('[%s]' % self.name, end='')
        print('Elapsed: %s' % (time.time() - self.tstart))

import json
class JsonDict(dict):
    """Provide json pretty-printing"""
    def __str__(self): return repr(self)
    def __repr__(self):
        s = json.dumps(self, indent=4, sort_keys=True, default=str)
        crop = lambda t: t[:80] + ("" if len(t)<80 else "...")
        s = "\n".join([crop(ln) for ln in s.split("\n")])
        return s
