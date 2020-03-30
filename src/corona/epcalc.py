##
"""Reproduce ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""

# TODO: evidences sum-to-1
# TODO: make SEIRS

## Imports
from corona.utils import *
from corona.model import *

## Model
model = SEIR2()

## Integrate

## Initial conditions, x0
x0 = zeros(10)
x0[1] = 0           # Exposed
x0[2] = 1/model.N   # Infected
x0[3] = 0           # Recovered
x0[0] = 1 - sum(x0) # Susceptible

# Time params
t_end = 200 # days
dt = 0.1
tt = linspace(0, t_end, int(t_end/dt)+1)

xx = zeros(tt.shape+x0.shape)
for k,t in enumerate(tt):
    if k: xx[k] = model.step(xx[k-1], t, t-tt[k-1])
    else: xx[k] = x0

# from scipy.integrate import odeint
# xx = odeint(model.dxdt, x0, tt)

# Multiply by population
xx = model.N * xx

# Facilitate unpacking
state = model.NamedState(*xx.T)

## Plot
fig, ax = freshfig(1)

# Normal plot:
# lbl='Susceptible'; ax.plot(tt, state.Susceptible, label=lbl, c=colrs[lbl])
# lbl='Exposed'    ; ax.plot(tt, state.Exposed, label=lbl, c=colrs[lbl])
# lbl='Infected'   ; ax.plot(tt, state.Infected, label=lbl, c=colrs[lbl])


# Barchart:
barchart = StackedBarChart(ax,state,tt)
barchart.add("Fatalities")
barchart.add("Hospitalized")
# barchart.add("Recovered")
barchart.add("Infected")
barchart.add("Exposed")
# barchart.add("Susceptible")


# Plot number of respirators
# nrk.no/vestland/mener-helsemyndighetene-overdriver-intensivkapasiteten-i-norge-1.14938514
nRespirators = 2000
ax.plot([0,t_end], 2*[nRespirators],"k--",lw=1,label="_nolegend_")
# ax.text(t_end, nRespirators,"Num. of respirators", va="bottom", ha="right",fontsize="small")
ax.text(0, nRespirators,"Num. of respirators", va="bottom", fontsize="small")

# Plot intervention line:
axY = ax.get_ylim()[1]
ax.plot(2*[model.intervention_time], [0,axY], "k--", lw=1,label="_nolegend_")
ax.text(model.intervention_time, axY,"Stricter measures", va="top", ha="right",fontsize="small",rotation=90)
