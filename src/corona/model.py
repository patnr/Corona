"""Model from ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""
from corona.utils import *
from collections import namedtuple


@dcs.dataclass
class SEIR2:
    # Params -- SEI (transmission dynamics)
    Rep          : float = 2.2 # Reproduction number
    # Timescales, where dt_X is the "mean" duration spent in state X.
    D_incbation  : float = 5.2 # Incubation dt
    D_infectious : float = 2.9 # Infection dt

    # Params -- IQR (clinical dynamics)
    # Proportions
    pDead             : float = 0.02  # Case fatality rate (CFR)
    pSevr             : float = 0.2   # Hospitalization
    # Note: Some of these timescales get compensated for previous states.
    Time_to_death     : float = 32    # Death
    D_recovery_mild   : float = 14    # Recovery mild
    D_recovery_severe : float = 31.5  # Recovery severe
    D_hospital_lag    : float = 5     # Time to hospitalization

    # Intervention params
    intervention_efficacy : float = 2/3 # slider parameter
    intervention_time     : float = 100 # day of stricter measures

    def __post_init__(self):

        # Aliases
        self.dt_I = self.D_infectious
        self.dt_E = self.D_incbation
        # Other parameterization: 
        # beta := Rep/dt_I : "Contact rate" (people/days)
        # a    :=   1/dt_E : "Incubation rate"
        # gamma:=   1/dt_I : "Recovery rate"

        # Corrections and aliases
        self.dt_mild = self.D_recovery_mild   - self.dt_I
        self.dt_fatl = self.Time_to_death     - self.dt_I
        self.dt_hosp = self.D_hospital_lag
        self.dt_sevr = self.D_recovery_severe - self.dt_I

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
        "Q_mild"     ,  # [3]
        "Q_sevr"     ,  # [4]
        "Q_hosp"     ,  # [5]
        "Q_fatl"     ,  # [6]
        # Recovered subgroups
        "R_mild"     ,  # [7]
        "R_sevr"     ,  # [8]
        "R_fatl"     ,  # [9]
        ])
    # Add diagnostic (non-prognostic) variables:
    class NamedVars(NamedState):
        @property
        def Hospitalized(self): return self.Q_hosp + self.Q_fatl
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

        # ------ Intervention switch ------
        if t>self.intervention_time:
            xRep = self.intervention_multiplier
        else:
            xRep = 1

        # ------ SEI (transmission dynamics) ------
        # Fluxes
        S2E = x.Susceptible * x.Infected / self.dt_I * self.Rep * xRep
        E2I = x.Exposed                  / self.dt_E
        I2Q = x.Infected                 / self.dt_I
        # Changes
        dS = -S2E
        dE = +S2E - E2I
        dI = -I2Q + E2I
        # dR = +I2Q

        # ------ IQR (clinical dynamics) ------
        # Fluxes
        Q2R_mild = x.Q_mild / self.dt_mild
        Q2R_sevr = x.Q_sevr / self.dt_hosp 
        Q2R_fatl = x.Q_fatl / self.dt_fatl        
        # Changes to Q
        dQ_mild = self.pMild*I2Q - Q2R_mild
        dQ_sevr = self.pSevr*I2Q - Q2R_sevr
        dQ_fatl = self.pDead*I2Q - Q2R_fatl
        # Changes to R
        dR_mild = +Q2R_mild
        dR_sevr = +Q2R_sevr
        dR_fatl = +Q2R_fatl

        # Hospitalized
        d_hosp = x.Q_sevr/self.dt_hosp - x.Q_hosp/self.dt_sevr # NB: wtf

        return np.asarray([dS, dE, dI, dQ_mild, dQ_sevr, d_hosp, dQ_fatl, dR_mild, dR_sevr, dR_fatl])
