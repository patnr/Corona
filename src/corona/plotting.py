"""Plotting tools.
"""
from corona.utils import *
import matplotlib as mpl
import matplotlib.pyplot as plt 
plt.ion()


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

thousands = mpl.ticker.StrMethodFormatter('{x:,.7g}')



leg_kws = dict(loc="upper left", bbox_to_anchor=(0.1,1), fontsize="8")

def reverse_legend(ax,**kws):
    "Reverse order of legend items in ``ax``."
    leg = ax.get_legend_handles_labels()
    leg = list(map(list, zip(*leg)))[::-1]
    ax.legend(*zip(*leg),**kws)



from matplotlib.widgets import Button, CheckButtons
def add_log_toggler(ax):
    def toggle_scale(_):
        # Get current status
        if isinstance(ax.toggle_log, mpl.widgets.Button):
            log_is_on = getattr(ax,"_log_is_on",False)
            ax._log_is_on = not log_is_on
        elif isinstance(ax.toggle_log, mpl.widgets.CheckButtons):
            log_is_on = not ax.toggle_log.get_status()[0]
        else: raise TypeError

        # Toggle
        if log_is_on:
            ax.set_yscale("linear")
            ax.yaxis.set_major_formatter(thousands)
        else:
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(thousands)
        plt.draw()
    # Toggle log-scale (button)
    toggle_ax = ax.figure.add_axes([0.8, 0.8, 0.15, 0.2], frameon=False)
    ax.toggle_log = CheckButtons(toggle_ax, ["Log scale"], [False])
    # toggle_ax = ax.figure.add_axes([0.8, 0.9, 0.15, 0.05])
    #ax.toggle_log = Button(toggle_ax, 'Toggle scale')
    ax.toggle_log.on_clicked(toggle_scale)



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
