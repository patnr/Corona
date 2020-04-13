"""Model from ``gabgoh.github.io/COVID``,

which uses a SEIR model elaborated with clinical dynamics.
"""
from corona.utils import *
from corona.maths import *


@dcs.dataclass
class SEIR2:
    # ------ Params SEI (transmission dynamics) ------
    Rep0            : float = 2.2 # Reproduction number
    # Intervention
    t_intervention  : float = 100 # Day of stricter measures
    dt_intervention : float = 30  # Duration of measures
    Rep_intervention: float = 0.9 # Reproduction number with stricter measures
    # Timescales, where dt_X is the "mean" duration spent in state X.
    dt_E            : float = 5.2 # Incubation (but not yet infect./sympt.)
    dt_I            : float = 2.9 # Infection (but not yet quaranteed)
    # Alternative parameterization: 
    # beta  := Rep0/dt_I : "Contact rate" (people/days)
    # a     :=    1/dt_E : "Incubation rate"
    # gamma :=    1/dt_I : "Recovery rate"

    # ------ Params IQR (clinical dynamics) ------
    # Proportions
    pDead   : float = 0.02  # Case fatality rate (CFR)
    pSevr   : float = 0.2   # Hospitalization
    @property               # Proportions must sum to 1, so:
    def pMild(self): return 1 - self.pSevr - self.pDead

    # Timescales, where dt_X is the "mean" duration spent in state X.
    dt_Q_mild : float = 14   - 2.9 # dt to recovery -- mild (since quarantine)
    dt_Q_fatl : float = 32   - 2.9 # dt to fatality         (since quarantine)
    dt_Q_sevr : float = 5          # dt to hospitalization  (since quarantine)
    dt_H      : float = 31.5 - 2.9 # dt to recovery (sevr)  (since hospitalization)
    # Note: gabgoh doesn't do dt_Q_sevr-=dt_I, nor dt_H-=(dt_I+dt_Q_sevr),
    #       which I think is wrong (from his definition of the timescales,
    #       which start from "showing symptoms").
    #       But his default param. values are left in place for reprodicubility.
    
    # ------ Variables ------
    # Use a namedtuple to unpack state variables.
    # Q: Why not just list them as args to dxdt?
    # A1: Make them accessible outside.
    # A2: dxdt(x,t) is neater than dxdt(a,b,c,...,t).
    # A3: dataclasses proved far too slow.
    Variables = namedtuple("PrognosticVars",[
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
    # Convenience:
    def init_state(self,**kwargs):
        """Init with kwargs and set susceptible = 1-everything_else."""
        x = {**{k:0 for k in self.Variables._fields}, **kwargs}
        x = np.array(list(x.values()))
        x[0] = 1 - sum(x[1:])
        return x

    def dxdt(self, state, t):
        "Dynamics."
        x = self.Variables(*state.T)

        # ------ Intervention switch ------
        if self.t_intervention < t < (self.t_intervention+self.dt_intervention):
            beta = self.Rep_intervention
        else:
            beta = self.Rep0

        # ------ SEI (transmission dynamics) ------
        # Fluxes
        S2E = x.Susceptible * x.Infected / self.dt_I * beta
        E2I = x.Exposed                  / self.dt_E
        I2Q = x.Infected                 / self.dt_I
        # Changes
        dS = -S2E
        dE = +S2E - E2I
        dI = -I2Q + E2I
        # dR = +I2Q

        # ------ IQR (clinical dynamics) ------
        # Fluxes
        Q2R_mild = x.Q_mild / self.dt_Q_mild
        Q2R_fatl = x.Q_fatl / self.dt_Q_fatl
        Q2H      = x.Q_sevr / self.dt_Q_sevr
        H2R      = x.Q_hosp / self.dt_H
        # Changes to Q
        dQ_mild = -Q2R_mild + self.pMild*I2Q 
        dQ_fatl = -Q2R_fatl + self.pDead*I2Q 
        dQ_sevr = -Q2H      + self.pSevr*I2Q 
        dH      = +Q2H      - H2R
        # Changes to R
        dR_mild = +Q2R_mild
        dR_fatl = +Q2R_fatl
        dR_sevr = +H2R

        return np.asarray([dS, dE, dI, dQ_mild, dQ_sevr, dH, dQ_fatl, dR_mild, dR_sevr, dR_fatl]).T


def with_diagnostics(cls):
    """Add some diagnostical (using the prognostic) variables."""
    class SubCls(cls):
        @property
        def Hospitalized(self): return self.Q_hosp + self.Q_fatl
        @property
        def Recovered(self): return self.R_mild + self.R_sevr
        @property
        def Fatalities(self): return self.R_fatl
    return SubCls
