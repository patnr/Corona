"""Plotting tools.
"""
from corona.utils import *
from corona.maths import *

import matplotlib as mpl
import matplotlib.pyplot as plt 
plt.ion()


from datetime import datetime, timedelta, timezone
import matplotlib.dates as mdates


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

# Coloschemes - https://www.schemecolor.com
cs_corona = "#6E1B09", "#D22C2C", "#F07249", "#5D5F5C", "#393A3C"
cs_aretro = "#EC5E64", "#6B3979", "#28A98F", "#FAD542", "#2ABBDA"
cs_pastel = "#998AD3", "#E494D3", "#CDF1AF", "#87DCC0", "#88BBE4"
cs_beach  = "#42B7C2", "#8FC8C4", "#FDF2C5", "#DECA98", "#A0795F", "#623D45"

colrs = dict(
        Fatalities="#386cb0",
        Hospitalized="#8da0cb",
        Recovered="#4daf4a",
        Infected="#f0027f",
        Exposed="#fdc086",
        Susceptible="grey",
        )

import hashlib
def colrz(component):
    if component in colrs:
        c = colrs[component]
    else:
        x = str(component).encode() # hashable
        # HASH = hash(tuple(x)) # Changes randomly each session
        HASH = int(hashlib.sha1(x).hexdigest(),16)
        colors = plt.get_cmap('tab20').colors
        c = colors[HASH%len(colors)]
    return c

thousands = mpl.ticker.StrMethodFormatter('{x:,.7g}')



leg_kws = dict(loc="upper left", bbox_to_anchor=(0.1,1), fontsize="8")

def reverse_legend(ax,**kws):
    "Reverse order of legend items in ``ax``."
    leg = ax.get_legend_handles_labels()
    leg = list(map(list, zip(*leg)))[::-1]
    ax.legend(*zip(*leg),**kws)



from matplotlib.widgets import Button, CheckButtons
def add_log_toggler(ax):
    """Add button that toggles log. scale."""

    def toggle_scale(_):
        "Toggle log. scale."
        # Get current status
        if isinstance(ax.toggle_log, mpl.widgets.Button):
            log_is_on = getattr(ax,"_log_is_on",False)
            ax._log_is_on = not log_is_on

            if log_is_on:
                # ax.toggle_log.label.set_text("Log. scale: Off")
                ax.toggle_log.color = "w"
            else:
                # ax.toggle_log.label.set_text("Log. scale: On")
                ax.toggle_log.color = "#D5F2E8"

        elif isinstance(ax.toggle_log, mpl.widgets.CheckButtons):
            log_is_on = not ax.toggle_log.get_status()[0]
        else: raise TypeError

        # Toggle
        if log_is_on:
            ax.set_yscale("linear")
            ax.set_ylim(bottom=0)
            ax.yaxis.set_major_formatter(thousands)
        else:
            ax.autoscale(True,axis="y")
            ax.set_yscale("log")
            ax.yaxis.set_major_formatter(thousands)
        plt.draw()

    # Get rectangle (position) for button placement,
    # but in figure coordinates, as required by fig.add_axes().
    height = .05
    width  = .15
    fig_coords = ax.figure.transFigure.inverted()

    # Place below legend
    # plt.pause(0.1) # Must draw before bbox of legend can be known.
    # leg = ax.get_legend()
    # bbox = leg.get_window_extent()
    # frame = leg.get_frame().get_bbox()
    # x,y,w,h = bbox.transformed(fig_coords).bounds
    # rect = [x, y-1.01*height, width, height]

    # Place in upper-right corner of ax
    x,y,w,h = ax.get_position().bounds
    rect = [x+w-width, y+h-1.01*height, width, height]

    # toggle_ax = ax.figure.add_axes(rect,frameon=True)
    # ax.toggle_log = Button(toggle_ax, 'Toggle scale',color="w")

    toggle_ax = ax.figure.add_axes(rect,frameon=False)
    ax.toggle_log = CheckButtons(toggle_ax, ["Log scale"], [False])
    dh = .3
    for box,cross in zip(ax.toggle_log.rectangles,ax.toggle_log.lines):
        box.set_y(dh)
        box.set_height(1-2*dh)
        cross[0].set_ydata([dh,1-dh])
        cross[1].set_ydata([dh,1-dh][::-1])


    ax.toggle_log.on_clicked(toggle_scale)


class CompartmentalPlot:
    def __init__(self,ax,state,tt,date0=None,**kwargs):
        self.ax        = ax
        self.state     = state
        self.tt        = tt
        self.kwargs    = kwargs

        if date0:
            if not date0.tzinfo:
                # MPL returns tz-aware (from user input, eg.),
                # so we need to make dates tz-aware for comparisons.
                date0 = date0.replace(tzinfo=timezone.utc)
            # Time-2-date conversion
            t2d = lambda t: date0 + (t-tt[0])*timedelta(1)
        else:
            # Passthrough
            t2d = lambda t: t
        self.date0 = date0
        self.t2d   = t2d
        self.dates = t2d(tt)

    def finalize(self):
        ax = self.ax

        # Adjust plot properties
        # ax.set_xlabel('Time (days)')
        # ax.set_ylabel('People')
        ax.set_ylim(bottom=0)
        ax.set_xlim(*self.dates[[0,-1]])

        # xticks -- datetime
        if self.date0:
            ax.xaxis.set_major_formatter(mdates.DateFormatter('%b'))
            ax.xaxis.set_major_locator(mdates.MonthLocator())
            ax.xaxis.set_minor_formatter(mdates.DateFormatter('%d'))
            # ax.xaxis.set_minor_locator(mdates.DayLocator(interval=10))
            ax.xaxis.set_minor_locator(mdates.AutoDateLocator())
            ax.figure.autofmt_xdate(bottom=0.11)

        # More adjustments:
        for edge in ["right","left","top"]:
            ax.spines[edge].set_visible(False)
        ax.grid(axis="y",ls=":",alpha=0.2, color="k")
        ax.yaxis.set_major_formatter(thousands)
        ax.tick_params(axis="y",pad=-1,length=0)
        ax.tick_params(axis="both",which="both",labelsize="small")
        ax.tick_params(axis="x",which="minor",labelsize="xx-small")
        plt.setp(ax.get_yticklabels(), alpha=0.4, ha="left", va="bottom", )
        plt.setp(ax.get_xticklabels(which="both"), alpha=0.4)
        # Would have to use axisartist for access to set_va, set_ha.

        add_log_toggler(ax)

        try:    __IPYTHON__
        except: plt.show(block=True)


class Lines(CompartmentalPlot):
    def add(self,label,**kwargs):
        yy = getattr(self.state, label)
        self.ax.plot(self.dates, yy.T, label=label, c=colrs[label],
                **{**self.kwargs, **kwargs})


class StackedBars(CompartmentalPlot):
    """A bar chart (histogram),

    but stacking each series on top of the previous.
    """
    def __init__(self,*args,**kwargs):
        super().__init__(*args,**kwargs)

        # Init stack
        self.stack = {}
        self.handles = {}

        # Plot resolution (in days)
        self.dt_plot = 2
        self.tt_plot = arange(*self.tt[[0,-1]],self.dt_plot)
        self.dates_plot = self.t2d(self.tt_plot)

        # Bar selection with highlighting and day info in legend
        self.alpha = .6
        self.ax.figure.canvas.mpl_connect('pick_event', self.onpick)

    def add(self,label):
        # Down-sample (interpolate)
        yy = np.interp(self.tt_plot, self.tt, getattr(self.state,label))
         
        # Accumulate bars
        cum = np.sum([y for y in self.stack.values()], 0)

        # Plot
        hh = self.ax.bar(self.dates_plot, yy, .8*self.dt_plot, bottom=cum,
                label=label, color=colrz(label),
                alpha=self.alpha, align="edge",picker=1)

        # Append bar heights to stack
        self.handles[label] = hh
        self.stack[label] = yy

        self.ax.legend(**leg_kws)
        reverse_legend(self.ax,**leg_kws)

    def onpick(self,event):
        rectangle = event.artist
        time = rectangle.xy[0]
        if self.date0:
            time = mdates.num2date(time)

        def arr_idx(arr,val):
            return abs(arr - val).argmin()

        def set_legend_for_day(t):
            iDay = arr_idx(self.dates, t)
            handles, labels = self.ax.get_legend_handles_labels()
            # Parse labels, insert numbers
            for i,lbl in enumerate(labels):
                num  = getattr(self.state,lbl)[iDay]
                new  = lbl.split(":")[0] + ": "
                new += thousands(round2sigfig(num,3))
                labels[i] = new
            # Get day (title) string
            day = int(self.tt[iDay])
            title = f"Day {day}"
            if self.date0:
                title += ": " + self.dates[iDay].strftime("%b %d, %Y")
            # Set legend
            self.ax.legend(handles[::-1],labels[::-1],
                    title=title, **leg_kws)
            plt.pause(0.2)

        def set_alpha_for_day(t):

            def setter(iDay,alpha):
                for label, rectangles in self.handles.items():
                    rectangles[iDay].set_alpha(alpha)

            # Reset alpha
            try:
                setter(self._iDay_alpha, self.alpha)
            except AttributeError:
                pass

            # Set alpha
            iDay = arr_idx(self.dates_plot, t)
            setter(iDay,1)
            plt.pause(0.01)
            self._iDay_alpha = iDay

        set_legend_for_day(time)
        set_alpha_for_day(time)
