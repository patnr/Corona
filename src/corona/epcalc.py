##
"""Reproduce ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.

The same code can also be used to run
the "simple" SEIR model, and also SIR.
"""

use_simple_model = False

# TODO: plot that evidences sum-to-1
# TODO: make Intervention
# TODO: make SEIRS

## Imports
from corona.utils import *

## Params
# Total population
N = 7*10**6

## Initial conditions, x0
x0 = zeros(4) if use_simple_model else zeros(10)
x0[1] = 0           # Exposed
x0[2] = 1/N         # Infected
x0[3] = 0           # Recovered
x0[0] = 1 - sum(x0) # Susceptible

## Models
if use_simple_model:
    @ens_compatible
    def dxdt(state, t):
        S,E,I,R = state
        dS = -beta*S*I
        dE = +beta*S*I - a*E
        dI = -gamma*I  + a*E
        dR = +gamma*I
        return dS, dE, dI, dR

    # Params
    # Note: set Tinc<<1 to turn SEIR into SIR.
    Rep  = 2 # Reproduction ratio, eg 2
    Tinc = 5 # Infection dt      , eg 5
    Tinf = 3 # Incubation dt     , eg 3
    # Proxy params:
    beta  = Rep/Tinf # Contact rate
    a     = 1/Tinc   # Incubation rate
    gamma = 1/Tinf   # Recovery rate

else:
    @ens_compatible
    def dxdt(state, t):
        S,E,I,  IM,IS,ISH,IF,  RM,RS,RF   =   state
        # N = sum(state)

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




## Integrate
t_end = 200 # days
dt = 0.1
tt = linspace(0, t_end, int(t_end/dt)+1)

# from scipy.integrate import odeint
# xx = odeint(dxdt, x0, tt)

xx = zeros(tt.shape+x0.shape)
for k,t in enumerate(tt):
    if k: xx[k] = rk4(dxdt, xx[k-1], t, t-tt[k-1])
    else: xx[k] = x0

# Unpack
SS, EE, II, *others = N * xx.T
if use_simple_model:     RR, = others
else: IM,IS,ISH,IF, RM,RS,RF = others

from dataclasses import dataclass
@dataclass
class NamedState:
    Susceptible : float
    Exposed     : float
    Infected    : float
    I_Mild      : float
    I_Sevr      : float
    I_Hosp      : float
    I_Fatl      : float
    R_Mild      : float
    R_Sevr      : float
    R_Fatl      : float
    @property
    def Fatalities(self): return self.R_Fatl
    @property
    def Hospitalized(self): return self.I_Hosp + self.I_Fatl

state = NamedState(*N*xx.T)

# Print #infected
get_k = lambda t: abs(tt-100).argmin()
sometime = 100
print(f"Infected at t={sometime} days:", int(round(II[get_k(sometime)])))



## Plot
plt.ion()
fig, ax = freshfig(1)

# Normal plot:
# # lbl='Susceptible'; ax.plot(tt, SS, label=lbl, c=colrs[lbl])
# lbl='Exposed'    ; ax.plot(tt, EE, label=lbl, c=colrs[lbl])
# lbl='Infected'   ; ax.plot(tt, II, label=lbl, c=colrs[lbl])
# if use_simple_model:
#     lbl='Recovered'; ax.plot(tt, RR, label=lbl, c=colrs[lbl])

def stack_bar(label):
    """Add a bar chart (histogram),

    but stack on top of previously added.
    """
    stack = stack_bar.stack
    ax    = stack_bar.ax

    # Down-sample (interpolate)
    _dt = 2 # plot resolution (in days)
    _xx = arange(0,t_end,_dt)
    _yy = np.interp(_xx, tt, getattr(state,label))
     
    cum = np.sum([y for y,l in stack], 0)

    h = ax.bar(_xx, _yy, .75*_dt, label=label, bottom=cum,
            color=colrs[label], alpha=0.7, align="edge",picker=5)
    stack_bar.handles += h

    stack.append((_yy,label))

# Init stacked bar chart
stack_bar.stack = []
stack_bar.ax = ax
stack_bar.handles = []
# Add bars
stack_bar("Fatalities")
stack_bar("Hospitalized")
stack_bar("Infected")
stack_bar("Exposed")


def todays_legend(day):
    handles, labels = ax.get_legend_handles_labels()
    for i,lbl in enumerate(labels):
        num  = getattr(state,lbl)[day]
        new  = lbl.split(":")[0]
        new += ": %d"%int(round2sigfig(num,3))
        labels[i] =new
    ax.legend(handles[::-1],labels[::-1], title="Day %d"%tt[day])
    plt.pause(0.01)


def onpick(event):
    rectangle = event.artist
    today = rectangle.xy[0]
    today_k = abs(tt-today).argmin()
    todays_legend(today_k)
fig.canvas.mpl_connect('pick_event', onpick)


# nrk.no/vestland/mener-helsemyndighetene-overdriver-intensivkapasiteten-i-norge-1.14938514
# nRespirators = 1000
# ax


# Adjust plot properties
ax.set_xlabel('Time (days)')
# ax.set_ylabel('People')
ax.legend()
ax.set_title("Click bars for number info.")
reverse_legend(ax)
# ax.set_ylim(0,9e5)
ax.set_xlim(0,t_end)
# Plot intervention line:
ax.plot(2*[InterventionTime], [0,ax.get_ylim()[1]], "k--", lw=1,label="_nolegend_")

# More adjustments:
for edge in ["right","left","top"]:
    ax.spines[edge].set_visible(False)
ax.grid(axis="y",ls="--",alpha=0.2, color="k")
thousands = mpl.ticker.StrMethodFormatter('{x:,.0f}')
ax.yaxis.set_major_formatter(thousands)
ax.tick_params(axis="y",pad=-1,length=0)
ax.tick_params(axis="both",labelsize="small")
_ = ax.get_yticklabels()
plt.pause(0.1) # avoid disappearing ticks bug
ax.set_yticklabels(_, ha="left", va="bottom")

##

##
