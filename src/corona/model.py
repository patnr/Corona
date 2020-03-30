"""Model from ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""
from corona.utils import *
from dataclasses import dataclass

@dataclass
class SEIR2:
    # Total population
    N : int = 7*10**6

    # Params -- Transmission dynamics
    Rep          : float = 2.2 # Reproduction number
    D_infectious : float = 2.9 # Infection dt (mean)
    D_incbation  : float = 5.2 # Incubation dt

    # Params -- Clinical dynamics
    # Durations:
    Time_to_death     : float = 32    # Death
    D_recovery_mild   : float = 14    # Recovery mild
    D_recovery_severe : float = 31.5  # Recovery severe
    D_hospital_lag    : float = 5     # Hospitalization
    # Proportions:
    CFR               : float = 0.02  # Case fatality
    pSevr             : float = 0.2   # Hospitalization

    # Intervention params
    intervention_efficacy : float = 2/3 # "to decrease transmission by"
    intervention_time     : float = 100 # day of stricter measures

    def __post_init__(self):
        # Aliases
        self.pDead   = self.CFR
        self.D_death = self.Time_to_death

        # Corrections
        self.D_recovery_mild   -= self.D_infectious
        self.D_recovery_severe -= self.D_infectious
        self.D_death           -= self.D_infectious
    @property
    def pMild(self): return 1 - self.pSevr - self.pDead
    @property
    def intervention_multiplier (self): return 1 - self.intervention_efficacy

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

    def step(self,x,t,dt):
        return rk4(self.dxdt, x, t, dt)

    # @ens_compatible
    def dxdt(self, state, t):

        # Proxy params
        beta  = self.Rep/(self.D_infectious) # Contact rate
        a     = 1/self.D_incbation           # Incubation rate
        gamma = 1/self.D_infectious          # Recovery rate (in people/days)

        S,E,I,  IM,IS,ISH,IF,  RM,RS,RF   =   state
        # IM  // Recovering (Mild)     
        # IS  // Recovering (Severe at home)
        # ISH // Recovering (Severe in hospital)
        # IF  // Recovering (Fatal)

        # RM // Recovered mild
        # RS // Recovered severe
        # RF // Fatal

        # Switch case
        b = beta
        if t>self.intervention_time:
            b *= self.intervention_multiplier

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
        dIR_Mild = IM / self.D_recovery_mild
        dIR_Sevr = IS / self.D_hospital_lag 
        dIR_Dead = IF / self.D_death        
        # Changes to Infected
        dMild   = self.pMild*I2R - dIR_Mild
        dSevr   = self.pSevr*I2R - dIR_Sevr
        dFatal  = self.pDead*I2R - dIR_Dead
        # Hospitalized
        dSevr_H = (1/self.D_hospital_lag)*IS - (1/self.D_recovery_severe)*ISH

        # RecoverED (or dead)
        dR_Mild   = +dIR_Mild
        dR_Sevr   = +dIR_Sevr
        dR_Fatal  = +dIR_Dead

        return np.asarray([dS, dE, dI, dMild, dSevr, dSevr_H, dFatal, dR_Mild, dR_Sevr, dR_Fatal])
