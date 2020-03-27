##
"""Reproduce gabgoh.github.io/COVID,

which uses a SEIR model elaborated with clinical dynamics.

The same code can also be used to run
the "simple" SEIR model, and also SIR.
"""

# TODO: recovery from death?
# TODO: plot that evidences sum-to-1
# TODO: make Intervention
# TODO: make SEIRS

## Imports
from corona.utils import *

## Params
use_simple_model = False
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
        # Fluxes between compartments
        dSE = b*S*I
        dEI = a*E
        dIR = gamma*I
        # Changes (for each compartment)
        dS = -dSE
        dE = +dSE - dEI
        dI = -dIR + dEI
        dR = +dIR # proxy var

        # Clinical dynamics
        # Fluxes to recovery (or death)
        dIR_Mild = IM / D_recovery_mild
        dIR_Sevr = IS / D_hospital_lag 
        dIR_Dead = IF / D_death        
        # Changes to Infected
        dMild   = pMild*dIR - dIR_Mild
        dSevr   = pSevr*dIR - dIR_Sevr
        dFatal  = pDead*dIR - dIR_Dead
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
dt = 1
tt = np.linspace(0, t_end, int(t_end/dt)+1)

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

# Print #infected
get_k = lambda t: abs(tt-100).argmin()
sometime = 100
print(f"Infected at t={sometime} days:", int(II[get_k(sometime)]))

## Plot
plt.ion()
fig, ax = freshfig(1)
# fatalities hostpitalized   Reocvered  Infectuous Exposed
# "#386cb0", "#8da0cb",      "#4daf4a", "#f0027f", "#fdc086"
# ax.plot(tt, SS, c="gray",    label='Susceptible')
ax.plot(tt, EE+II, c="#fdc086", label='Exposed (+I)')
ax.plot(tt, II, c="#f0027f", label='Infected')
if use_simple_model:
    ax.plot(tt, RR, c="#4daf4a", label='Recovered')
ax.set_xlabel('Time (days)')
ax.set_ylabel('People')
# ax.set_ylim(0,9e5)
ax.grid()
legend = ax.legend()

##

##
