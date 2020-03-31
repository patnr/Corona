##
"""Reproduce ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""

# TODO: make SEIRS

## Imports
from corona.utils import *
from corona.model import *
from corona.plotting import *

## Params

# model = SEIR2(Rep=3)
model = SEIR2()

# Population size
N = 7*10**6

# Time -- unit: days
t_end = 200
dt = 0.1
tt = linspace(0, t_end, int(t_end/dt)+1)

## Integrate
x0 = model.init_state(Infected=1/N)
xx = zeros(tt.shape+x0.shape)
for k,t in enumerate(tt):
    if k: xx[k] = model.step(xx[k-1], t, t-tt[k-1])
    else: xx[k] = x0

# from scipy.integrate import odeint
# xx = odeint(model.dxdt, x0, tt)

# Multiply by population
xx = N * xx


## Plot
fig, ax = freshfig(1)

state = model.NamedVars(*xx.T)

# Normal plot:
# # lbl='Susceptible'; ax.plot(tt, getattr(state,lbl), label=lbl, c=colrs[lbl])
# lbl='Exposed'    ; ax.plot(tt, getattr(state,lbl), label=lbl, c=colrs[lbl])
# lbl='Infected'   ; ax.plot(tt, getattr(state,lbl), label=lbl, c=colrs[lbl])


# Barchart:
barchart = StackedBarChart(ax,state,tt)
barchart.add("Hospitalized")
barchart.add("Fatalities")
# barchart.add("Recovered")
barchart.add("Infected")
barchart.add("Exposed")
# barchart.add("Susceptible")

# # All state variables (Note: sums to 1):
# for f in state._fields:
#     if f != "Q_hosp":
#         barchart.add(f)



add_log_toggler(ax)

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
