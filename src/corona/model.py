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




    # Use namedtuple to unpack state variables.
    # Q: Why not just list them as args to dxdt?
    # A: Make them accessible outside.
    # Note: using dataclasses was far too slow.
    NamedState = namedtuple("PrognosticVars",[
        "Susceptible",  # [0]
        "Exposed"    ,  # [1]
        "Infected"   ,  # [2]
        # Infected subgroups
        "I_mild"     ,  # [3]
        "I_sevr"     ,  # [4]
        "I_hosp"     ,  # [5]
        "I_fatl"     ,  # [6]
        # Recovered subgroups
        "R_mild"     ,  # [7]
        "R_sevr"     ,  # [8]
        "R_fatl"     ,  # [9]
        ])
    # Add diagnostic (non-prognostic) variables:
    class NamedVars(NamedState):
        @property
        def Hospitalized(self): return self.I_hosp + self.I_fatl
        @property
        def Recovered(self): return self.R_mild + self.R_sevr
        @property
        def Fatalities(self): return self.R_fatl
    # Convenience:
    def init_state(self,**kwargs):
        """Init with kwargs and set susceptible = 1-everything_else."""
        x = {**{k:0 for k in self.NamedState._fields}, **kwargs}
        x = np.array(list(x.values()))
        x[0] = 1 - sum(x[1:])
        return x




    def step(self,x,t,dt):
        "Integrate dxdt over dt."
        return rk4(self.dxdt, x, t, dt)

    # @ens_compatible
    def dxdt(self, state, t):
        "Dynamics."
        x = self.NamedState(*state)

        # Proxy params
        beta  = self.Rep/(self.D_infectious) # Contact rate
        a     = 1/self.D_incbation           # Incubation rate
        gamma = 1/self.D_infectious          # Recovery rate (in people/days)

        # Switch case
        b = beta
        if t>self.intervention_time:
            b *= self.intervention_multiplier

        # SEIR
        # Fluxes
        S2E = b*x.Susceptible*x.Infected
        E2I = a*x.Exposed
        I2R = gamma*x.Infected
        # Changes (for each compartment)
        dS = -S2E
        dE = +S2E - E2I
        dI = -I2R + E2I
        dR = +I2R # proxy var

        # Clinical dynamics
        # Fluxes to recovery (or death)
        dIR_mild = x.I_mild / self.D_recovery_mild
        dIR_sevr = x.I_sevr / self.D_hospital_lag 
        dIR_dead = x.I_fatl / self.D_death        
        # Changes to Infected
        d_mild = self.pMild*I2R - dIR_mild
        d_sevr = self.pSevr*I2R - dIR_sevr
        d_fatl = self.pDead*I2R - dIR_dead
        # Hospitalized
        d_hosp = (1/self.D_hospital_lag)*x.I_sevr - (1/self.D_recovery_severe)*x.I_hosp

        # RecoverED (or dead)
        dR_mild = +dIR_mild
        dR_sevr = +dIR_sevr
        dR_fatl = +dIR_dead

        return np.asarray([dS, dE, dI, d_mild, d_sevr, d_hosp, d_fatl, dR_mild, dR_sevr, dR_fatl])
