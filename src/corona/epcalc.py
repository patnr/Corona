##
"""Reproduce ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.

The same code can also be used to run
the "simple" SEIR model, and also SIR.
"""

# TODO: evidences sum-to-1
# TODO: make SEIRS

## Imports
from corona.utils import *

## Params

# Total population
N = 7*10**6

# Time params
t_end = 200 # days
dt = 0.1
tt = linspace(0, t_end, int(t_end/dt)+1)

# Params -- Transmission dynamics
Rep          = 2.2 # Reproduction number
D_infectious = 2.9 # Infection dt (mean)
D_incbation  = 5.2 # Incubation dt
# Proxies
beta  = Rep/(D_infectious) # Contact rate
a     = 1/D_incbation      # Incubation rate
gamma = 1/D_infectious     # Recovery rate (in people/days)

# Params -- Clinical dynamics
Time_to_death     = 32
D_recovery_mild   = 14            - D_infectious # Recovery dt mild
D_recovery_severe = 31.5          - D_infectious # Recovery dt severe
D_hospital_lag    = 5                            # Hospitalization dt
D_death           = Time_to_death - D_infectious # Death dt
CFR               = 0.02                         # Case fatality proportion
pDead             = CFR                          # Alias
pSevr             = 0.2                          # Hospitalization proportion
pMild             = 1 - pSevr - pDead            # Mild            proportion

# Intervention params
intervention_efficacy   = 2/3 # "to decrease transmission by"
intervention_multiplier = 1 - intervention_efficacy
InterventionTime        = 100




## Initial conditions, x0
x0 = zeros(10)
x0[1] = 0           # Exposed
x0[2] = 1/N         # Infected
x0[3] = 0           # Recovered
x0[0] = 1 - sum(x0) # Susceptible

## Model
@ens_compatible
def dxdt(state, t):
    S,E,I,  IM,IS,ISH,IF,  RM,RS,RF   =   state
    # IM  // Recovering (Mild)     
    # IS  // Recovering (Severe at home)
    # ISH // Recovering (Severe in hospital)
    # IF  // Recovering (Fatal)

    # RM // Recovered mild
    # RS // Recovered severe
    # RF // Fatal

    b = beta
    if t>InterventionTime:
        b *= intervention_multiplier

    # SEIR
    # Fluxes
    S2E = b*S*I
    E2I = a*E
    I2R = gamma*I
    # Changes (for each compartment)
    dS = -S2E
    dE = +S2E - E2I
    dI = -I2R + E2I
    dR = +I2R # proxy var

    # Clinical dynamics
    # Fluxes to recovery (or death)
    dIR_Mild = IM / D_recovery_mild
    dIR_Sevr = IS / D_hospital_lag 
    dIR_Dead = IF / D_death        
    # Changes to Infected
    dMild   = pMild*I2R - dIR_Mild
    dSevr   = pSevr*I2R - dIR_Sevr
    dFatal  = pDead*I2R - dIR_Dead
    # Hospitalized
    dSevr_H = (1/D_hospital_lag)*IS - (1/D_recovery_severe)*ISH

    # RecoverED (or dead)
    dR_Mild   = +dIR_Mild
    dR_Sevr   = +dIR_Sevr
    dR_Fatal  = +dIR_Dead

    return np.asarray([dS, dE, dI, dMild, dSevr, dSevr_H, dFatal, dR_Mild, dR_Sevr, dR_Fatal])



## Integrate
# from scipy.integrate import odeint
# xx = odeint(dxdt, x0, tt)

xx = zeros(tt.shape+x0.shape)
for k,t in enumerate(tt):
    if k: xx[k] = rk4(dxdt, xx[k-1], t, t-tt[k-1])
    else: xx[k] = x0
# Multiply by population
xx = N*xx



# Facilitate unpacking
from dataclasses import dataclass
@dataclass
class NamedState:
    Susceptible : float = 0
    Exposed     : float = 0 
    Infected    : float = 0 
    I_Mild      : float = 0 
    I_Sevr      : float = 0 
    I_Hosp      : float = 0 
    I_Fatl      : float = 0 
    R_Mild      : float = 0 
    R_Sevr      : float = 0 
    R_Fatl      : float = 0  
    # Non-prognostic variables:
    @property
    def Hospitalized(self): return self.I_Hosp + self.I_Fatl
    @property
    def Recovered(self): return self.R_Mild + self.R_Sevr
    @property
    def Fatalities(self): return self.R_Fatl
state = NamedState(*xx.T)

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
fig.canvas.mpl_connect('pick_event', barchart.onpick)

# Plot number of respirators
# nrk.no/vestland/mener-helsemyndighetene-overdriver-intensivkapasiteten-i-norge-1.14938514
nRespirators = 10000
ax.plot([0,t_end], 2*[nRespirators],"k--",lw=1,label="_nolegend_")
# ax.text(t_end, nRespirators,"Num. of respirators", va="bottom", ha="right",fontsize="small")
ax.text(0, nRespirators,"Num. of respirators", va="bottom", fontsize="small")

# Plot intervention line:
axY = ax.get_ylim()[1]
ax.plot(2*[InterventionTime], [0,axY], "k--", lw=1,label="_nolegend_")
ax.text(InterventionTime, axY,"Stricter measures", va="top", ha="right",fontsize="small",rotation=90)

# Adjust plot properties
ax.set_xlabel('Time (days)')
# ax.set_ylabel('People')
ax.legend()
reverse_legend(ax,**leg_kws)
# ax.set_ylim(0,9e5)
ax.set_xlim(0,t_end)

# More adjustments:
for edge in ["right","left","top"]:
    ax.spines[edge].set_visible(False)
ax.grid(axis="y",ls="--",alpha=0.2, color="k")
ax.yaxis.set_major_formatter(thousands)
ax.tick_params(axis="y",pad=-1,length=0)
ax.tick_params(axis="both",labelsize="small")
_ = ax.get_yticklabels()
plt.pause(0.1) # avoid disappearing ticks bug
ax.set_yticklabels(_, ha="left", va="bottom")

try:    __IPYTHON__
except: plt.show(block=True)

##

##
