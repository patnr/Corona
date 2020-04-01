##
"""Reproduce ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""

## Imports
from corona.utils import *
from corona.maths import *
from corona.model import *
from corona.plotting import *

## Params

# model = SEIR2(Rep=3)
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
state = model.NamedVars(*xx.T)


## Plot
fig, ax = freshfig(1)

# Normal plot:
# cPlot = Lines(ax,state,tt,date0)
# cPlot.add("Exposed")
# cPlot.add("Infected")

# Barchart:
cPlot = StackedBars(ax,state,tt,date0)
cPlot.add("Hospitalized")
cPlot.add("Fatalities")
# cPlot.add("Recovered")
cPlot.add("Infected")
cPlot.add("Exposed")
# cPlot.add("Susceptible")

# All state variables (should yield constant topline):
# for f in state._fields:
#     cPlot.add(f)

cPlot.finalize()

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
ax.plot(2*[cPlot.t2d(model.t_intervention)], yL, "k--", lw=1,label="_nolegend_")
ax.text(cPlot.t2d(model.t_intervention), yL[1],"Stricter measures", va="top", ha="right",fontsize="small",rotation=90)

##
##
