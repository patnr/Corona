"""Model from ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""
from corona.utils import *
from collections import namedtuple


@dcs.dataclass
class SEIR2:
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




    @dcs.dataclass
    class NamedState:
        Susceptible : float = 0
        Exposed     : float = 0 
        Infected    : float = 0 
        I_mild      : float = 0 
        I_sevr      : float = 0 
        I_hosp      : float = 0 
        I_fatl      : float = 0 
        R_mild      : float = 0 
        R_sevr      : float = 0 
        R_fatl      : float = 0  
        # Non-prognostic variables:
        @property
        def Hospitalized(self): return self.I_hosp + self.I_fatl
        @property
        def Recovered(self): return self.R_mild + self.R_sevr
        @property
        def Fatalities(self): return self.R_fatl
        # Convenience:
        def __post_init__(self):
            complement = [v for k,v in dcs.asdict(self).items() if k!="Susceptible"]
            self.Susceptible = 1 - sum(complement)
        def asarray(self):
            return np.asarray(dcs.astuple(self))






    def step(self,x,t,dt):
        return rk4(self.dxdt, x, t, dt)

    # @ens_compatible
    def dxdt(self, state, t):
        x = self.NamedState(*state)

        # Proxy params
        beta  = self.Rep/(self.D_infectious) # Contact rate
        a     = 1/self.D_incbation           # Incubation rate
        gamma = 1/self.D_infectious          # Recovery rate (in people/days)

        Susceptible,Exposed,Infected,  I_mild,I_sevr,I_hosp,I_fatl,  R_mild,R_sevr,R_fatl   =   state

        # Switch case
        b = beta
        if t>self.intervention_time:
            b *= self.intervention_multiplier

        # SEIR
        # Fluxes
        S2E = b*x.Susceptible*Infected
        E2I = a*Exposed
        I2R = gamma*Infected
        # Changes (for each compartment)
        dS = -S2E
        dE = +S2E - E2I
        dI = -I2R + E2I
        dR = +I2R # proxy var

        # Clinical dynamics
        # Fluxes to recovery (or death)
        dIR_Mild = I_mild / self.D_recovery_mild
        dIR_Sevr = I_sevr / self.D_hospital_lag 
        dIR_Dead = I_fatl / self.D_death        
        # Changes to Infected
        dMild   = self.pMild*I2R - dIR_Mild
        dSevr   = self.pSevr*I2R - dIR_Sevr
        dFatal  = self.pDead*I2R - dIR_Dead
        # Hospitalized
        dSevr_H = (1/self.D_hospital_lag)*I_sevr - (1/self.D_recovery_severe)*I_hosp

        # RecoverED (or dead)
        dR_Mild   = +dIR_Mild
        dR_Sevr   = +dIR_Sevr
        dR_Fatal  = +dIR_Dead

        return np.asarray([dS, dE, dI, dMild, dSevr, dSevr_H, dFatal, dR_Mild, dR_Sevr, dR_Fatal])
