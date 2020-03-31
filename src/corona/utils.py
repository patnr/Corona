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
plt.ion()

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


colrs = dict(
        Fatalities="#386cb0",
        Hospitalized="#8da0cb",
        Recovered="#4daf4a",
        Infected="#f0027f",
        Exposed="#fdc086",
        Susceptible="grey",
        )

thousands = mpl.ticker.StrMethodFormatter('{x:,.0f}')

leg_kws = dict(loc="upper left", bbox_to_anchor=(0.1,1), fontsize="8")

def reverse_legend(ax,**kws):
    "Reverse order of legend items in ``ax``."
    leg = ax.get_legend_handles_labels()
    leg = list(map(list, zip(*leg)))[::-1]
    ax.legend(*zip(*leg),**kws)


from matplotlib.widgets import Button, CheckButtons

class StackedBarChart:
    """A bar chart (histogram),

    but stacking each series on top of the previous.
    """
    def __init__(self,ax,state,tt_full):
        self.ax = ax
        self.state = state
        self.tt_full = tt_full

        self.stack = {}
        self.handles = {}

        self.dt = 2 # plot resolution (in days)
        t_end = tt_full[-1]
        self.tt = arange(0,t_end,self.dt)

        # Toggle log-scale (button)
        toggle_ax = ax.figure.add_axes([0.8, 0.8, 0.15, 0.2], frameon=False)
        self.toggle_log = CheckButtons(toggle_ax, ["Log scale"], [False])
        # toggle_ax = ax.figure.add_axes([0.8, 0.9, 0.15, 0.05])
        # self.toggle_log = Button(toggle_ax, 'Toggle scale')
        self.toggle_log.on_clicked(self.toggle_scale)

        # Highlight bars and add day info in legend
        self.alpha = .65
        ax.figure.canvas.mpl_connect('pick_event', self.onpick)

        # Adjust plot properties
        ax.set_xlabel('Time (days)')
        # ax.set_ylabel('People')
        # ax.set_ylim(0,9e5)
        ax.set_xlim(0,t_end)

        # More adjustments:
        for edge in ["right","left","top"]:
            ax.spines[edge].set_visible(False)
        ax.grid(axis="y",ls="--",alpha=0.2, color="k")
        ax.yaxis.set_major_formatter(thousands)
        ax.tick_params(axis="y",pad=-1,length=0)
        ax.tick_params(axis="both",labelsize="small")
        plt.setp(ax.get_yticklabels(), ha="left", va="bottom")
        # Would have to use axisartist for access to set_va, set_ha.

        try:    __IPYTHON__
        except: plt.show(block=True)

    def add(self,label):
        # Down-sample (interpolate)
        yy = np.interp(self.tt, self.tt_full, getattr(self.state,label))
         
        # Accumulate bars
        cum = np.sum([y for y in self.stack.values()], 0)

        # Plot
        hh = self.ax.bar(self.tt, yy, .6*self.dt, bottom=cum,
                label=label, color=colrs[label],
                alpha=self.alpha, align="edge",picker=5)

        # Append bar heights to stack
        self.handles[label] = hh
        self.stack[label] = yy

        self.ax.legend(**leg_kws)
        reverse_legend(self.ax,**leg_kws)

    def day_index(self,t):      return abs(self.tt      - t).argmin()
    def day_index_full(self,t): return abs(self.tt_full - t).argmin()

    def set_legend_for_day(self,t):
        iDay = self.day_index_full(t)
        handles, labels = self.ax.get_legend_handles_labels()
        for i,lbl in enumerate(labels):
            num  = getattr(self.state,lbl)[iDay]
            new  = lbl.split(":")[0] + ": "
            new += thousands(round2sigfig(num,3))
            labels[i] = new
        self.ax.legend(handles[::-1],labels[::-1],
                title="Day %d"%self.tt_full[iDay],
                **leg_kws)
        plt.pause(0.01)

    def set_alpha_for_day(self,t):

        def setter(iDay,alpha):
            for label, rectangles in self.handles.items():
                rectangles[iDay].set_alpha(alpha)

        # Reset alpha
        try:
            setter(self._iDay_alpha, self.alpha)
        except AttributeError:
            pass

        # Set alpha
        iDay = self.day_index(t)
        setter(iDay,1)
        plt.pause(0.01)
        self._iDay_alpha = iDay

    def onpick(self,event):
        rectangle = event.artist
        time = rectangle.xy[0]
        self.set_legend_for_day(time)
        self.set_alpha_for_day(time)

    def toggle_scale(self,_):
        # Get current status
        if isinstance(self.toggle_log, mpl.widgets.Button):
            log_is_on = getattr(self,"_log_is_on",False)
            self._log_is_on = not log_is_on
        elif isinstance(self.toggle_log, mpl.widgets.CheckButtons):
            log_is_on = not self.toggle_log.get_status()[0]
        else: raise TypeError

        # Toggle
        if log_is_on:
            self.ax.set_yscale("linear")
            self.ax.yaxis.set_major_formatter(thousands)
        else:
            self.ax.set_yscale("log")
            f = mpl.ticker.FuncFormatter(lambda y, _: '{:g}'.format(y))
            self.ax.yaxis.set_major_formatter(f)
        plt.draw()





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





def round2sigfig(num,nfig=1):
    """Round number to significant figures"""

    def ndecimal(x):
        if x==0 or not np.isfinite(x):
            # "Behaviour not defined" => should not be relied upon.
            return 1
        else:
            return -int(floor(log10(abs(x))))

    nfig =nfig-1
    n    = nfig + ndecimal(num)
    return np.round(num, n) # n specified => float (always)



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
