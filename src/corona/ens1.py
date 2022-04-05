"""
"""

## Imports
from corona.utils import *
from corona.maths import *
from corona.model import *
from corona.plotting import *

## Params
model = SEIR2(t_intervention=100)
# Population size
nPop = 5*10**6
# Ens size
N = 20
# Time -- unit: days
t_end = 365
dt    = 1
tt    = linspace(0 , t_end , int(t_end/dt)+1)
date0 = datetime(2020,2,12)
# date0 = None

## Init ensemble
x0 = model.init_state(Infected=1/nPop)
E0 = x0 + zeros((N,len(x0)))
E0[:,model.Variables._fields.index("Infected")] = .1/nPop * randn(N)
E0 = E0.clip(min=1e-9)

## Integrate
EE = nPop*integrate(model.dxdt, E0, tt)

# Unpack
NamedVars = with_diagnostics(model.Variables)
state = NamedVars(*EE.T)

## Plot
fig, ax = freshfig(1)

coPlot = Lines(ax,state,tt, date0)
coPlot.add("Exposed")
coPlot.add("Infected")
coPlot.add("Hospitalized")
coPlot.add("Fatalities")
coPlot.add("Recovered")
coPlot.finalize()
