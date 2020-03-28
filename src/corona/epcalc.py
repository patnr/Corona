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
    def Hospitalized(self): return self.I_Hosp + self.I_Fatl
    @property
    def Recovered(self): return self.R_Mild + self.R_Sevr
    @property
    def Fatalities(self): return self.R_Fatl

state = NamedState(*N*xx.T)




## Plot
plt.ion()
fig, ax = freshfig(1)

# Normal plot:
# # lbl='Susceptible'; ax.plot(tt, SS, label=lbl, c=colrs[lbl])
# lbl='Exposed'    ; ax.plot(tt, EE, label=lbl, c=colrs[lbl])
# lbl='Infected'   ; ax.plot(tt, II, label=lbl, c=colrs[lbl])
# if use_simple_model:
#     lbl='Recovered'; ax.plot(tt, RR, label=lbl, c=colrs[lbl])

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
        self.tt = arange(0,t_end,self.dt)

        self.alpha = .65

    def add(self,label):
        # Down-sample (interpolate)
        yy = np.interp(self.tt, self.tt_full, getattr(self.state,label))
         
        # Accumulate bars
        cum = np.sum([y for y in self.stack.values()], 0)

        # Plot
        hh = self.ax.bar(self.tt, yy, .75*self.dt, bottom=cum,
                label=label, color=colrs[label],
                alpha=self.alpha, align="edge",picker=5)

        # Append bar heights to stack
        self.handles[label] = hh
        self.stack[label] = yy

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



# Add bars
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
leg_kws = dict(loc="upper left", bbox_to_anchor=(0.1,1), fontsize="8")
ax.legend(
        list(barchart.handles.values())[::-1],
        list(barchart.handles.keys())[::-1],
        **leg_kws)
# reverse_legend(ax)
# ax.legend(**leg_kws)
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

##

##
