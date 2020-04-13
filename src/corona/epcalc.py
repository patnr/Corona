##
"""Reproduce ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""

## Imports
from corona.utils import *
from corona.maths import *
from corona.model import *
from corona.plotting import *
from mpl_tools import freshfig

## Params

# model = SEIR2(Rep0=3)
model = SEIR2()

# Population size
nPop = 7*10**6

# Time -- unit: days
t_end = 200
dt    = 0.1
tt    = linspace(0 , t_end , int(t_end/dt)+1)
date0 = datetime(2020,2,12)

## Integrate
x0 = model.init_state(Infected=1/nPop)
# xx = scipy.integrate.odeint(model.dxdt, x0, tt)
xx = integrate(model.dxdt, x0, tt)

# Multiply by population
xx = nPop * xx

# Unpack
NamedVars = with_diagnostics(model.Variables)
state = NamedVars(*xx.T)


## Plot
fig, ax = freshfig(1)

# Normal plot:
# coPlot = Lines(ax,state,tt,date0)
# coPlot.add("Exposed")
# coPlot.add("Infected")

# Barchart:
coPlot = StackedBars(ax,state,tt,date0)
coPlot.add("Hospitalized")
coPlot.add("Fatalities")
# coPlot.add("Recovered")
coPlot.add("Infected")
coPlot.add("Exposed")
# coPlot.add("Susceptible")

# All state variables (should yield constant topline):
# for f in state._fields:
#     coPlot.add(f)

coPlot.finalize()
reverse_legend(ax,**leg_kws)

## Add more info to plot
# Plot number of respirators
# nrk.no/vestland/mener-helsemyndighetene-overdriver-intensivkapasiteten-i-norge-1.14938514
nRespirators = 2000
xL = ax.get_xlim()
ax.plot(xL, 2*[nRespirators],"k--",lw=1,label="_nolegend_")
# ax.text(t_end, nRespirators,"Num. of respirators", va="bottom", ha="right",fontsize="small")
ax.text(xL[0], nRespirators,"Num. of respirators", va="bottom", fontsize="small")

# Plot intervention line:
yL = ax.get_ylim()
ax.plot(2*[coPlot.t2d(model.t_intervention)], yL, "k--", lw=1,label="_nolegend_")
ax.text(coPlot.t2d(model.t_intervention), yL[1],"Stricter measures", va="top", ha="right",fontsize="small",rotation=90)

##
##
