"""Misc
"""

# Already in python_startup.py
import sys
import os

# Others
import time
# import re
# import builtins
# import dataclasses as dcs
# from typing import Optional, Any


import numpy as np
import scipy as sp
import numpy.random
import scipy.linalg as sla
import numpy.linalg as nla
import scipy.stats as ss

from scipy.linalg import svd
from numpy.linalg import eig
# eig() of scipy.linalg necessitates using np.real_if_close().
from scipy.linalg import sqrtm, inv, eigh

from numpy import \
    pi, nan, \
    log, log10, exp, sin, cos, tan, \
    sqrt, floor, ceil, \
    mean, prod, \
    diff, cumsum, \
    array, asarray, asmatrix, \
    linspace, arange, reshape, \
    eye, zeros, ones, diag, trace \
    # Don't shadow builtins: sum, max, abs, round, pow

import matplotlib as mpl
import matplotlib.pyplot as plt 
# plt.ion()

np.set_printoptions(suppress=True,threshold=200,precision=6)
# Instead of set_np_linewidth, just let terminal do wrapping:
np.set_printoptions(linewidth=9999)




def freshfig(num=None,figsize=None,*args,**kwargs):
    """Create/clear figure.

    Similar to::

      fig, ax = suplots(*args,**kwargs)

    With the modification that:

    - If the figure does not exist: create it.
      This allows for figure sizing -- even with mpl backend MacOS.
    - Otherwise: clear figure.
      Avoids closing/opening so as to keep pos and size.
    """
    exists = plt.fignum_exists(num)

    fig = plt.figure(num=num,figsize=figsize)
    fig.clear()

    _, ax = plt.subplots(num=fig.number,*args,**kwargs)
    return fig, ax







def rk4(f, x, t, dt, order=4):
    """Runge-Kutta N-th order (explicit, non-adaptive) numerical ODE solvers.""" 
    if order >=1: k1 = dt * f( x        , t      )
    if order >=2: k2 = dt * f( x+k1/2   , t+dt/2 )
    if order ==3: k3 = dt * f( x+k2*2-k1, t+dt   )
    if order ==4:                         
                  k3 = dt * f( x+k2/2   , t+dt/2 )
                  k4 = dt * f( x+k3     , t+dt   )
    if    order ==1: return x + k1
    elif  order ==2: return x + k2
    elif  order ==3: return x + (k1 + 4*k2 + k3)/6
    elif  order ==4: return x + (k1 + 2*(k2 + k3) + k4)/6
    else: raise NotImplementedError


import functools
def ens_compatible(func):
    """Tranpose before and after.

    Helpful to make functions compatible with both 1d and 2d ndarrays.

    An older version also used np.atleast_2d and squeeze(),
    but that is more messy than necessary.

    Note: this is not the_way™ -- other tricks are sometimes more practical.
    See for example core.py:dxdt() of LorenzUV, Lorenz96, LotkaVolterra.
    """
    @functools.wraps(func)
    def wrapr(x,*args,**kwargs):
        return np.asarray(func(x.T,*args,**kwargs)).T
    return wrapr
