"""
"""

## Imports
from corona.utils import *
from corona.maths import *
from corona.model import *
from corona.plotting import *

np.random.seed(3)

## Time -- unit: days
t_end = 365
dt    = 1
tt    = linspace(0 , t_end , int(t_end/dt)+1)
date0 = datetime(2020,2,26) # 1st confirmed case in Norway
# date0 = None


## Params
# model = SEIR2(dt_E=5.2,dt_I=2.9,dt_Q_mild=14-2.9,dt_Q_sevr=5,dt_H=31.5-2.9,dt_Q_fatl=32,
        # pDead=0.01,pSevr=0.02,
        # Rep0=3.8,Rep_intervention=0.9,
        # t_intervention=15,dt_intervention=30)

model        = SEIR2(t_intervention=15,dt_intervention=30)
variables    = model.Variables._fields
parameters   = tuple(k for k in vars(model) if "t_inter" not in k)

StateVector = namedtuple("Estimated", variables + parameters)
NamedVars = with_diagnostics(StateVector)

# Population size
nPop = 5.368e6


## Init ensemble
def iChi2(mean,var=None):
    if var is None:
        var = (mean/10)**2
    def _iChi2(s,nu):
        return ss.invgamma(a=nu/2, scale=nu/2*s)
    nu = 4 + 2*mean**2/var
    s = mean*(nu-2)/nu
    return _iChi2(s,nu)

def X0_sample(N):
    # Defaults
    ens_var = {k:               0*ones(N) for k in variables}
    ens_par = {k:getattr(model,k)*ones(N) for k in parameters}
    E = {**ens_var, **ens_par}

    # E["Infected"] =  1/nPop + 0.1/nPop * randn(N)
    # E["Exposed"]  = 10/nPop +   1/nPop * randn(N)
    _mean = array([1,10]) / nPop
    _corr = array([[1,.8],[.8,1]])
    _std  = _mean / 10
    _cov  = _std[:,None] * _corr * _std
    IE = ss.multivariate_normal(_mean, _cov).rvs(N).clip(min=0)
    E["Infected"], E["Exposed"] = IE.T
    E["Susceptible"] = 1 - E["Infected"] - E["Exposed"]

    E["Rep0"             ] = iChi2(model.Rep0, 1.5**2).rvs(N)
    E["Rep_intervention" ] = iChi2(model.Rep_intervention, 0.7**2).rvs(N)
    E["pDead"            ] = iChi2(model.pDead).rvs(N)
    E["pSevr"            ] = iChi2(model.pSevr).rvs(N)

    assert np.all( E["pDead"] + E["pSevr"] < 1 )
    for k in E: assert np.all( 0 <= E[k] )
    return E

# Ens size
N = 300
E0 = X0_sample(N)
# sns.jointplot("pDead", "pSevr", data=E, marginal_kws=dict(bins=50, rug=True), kind="reg")
E = np.array(list(E0.values())).T


## Obs
# Load
obs_df = pd.read_csv("Norway.txt", sep='\s+', index_col=0, parse_dates=True)
# Compute days relative to date0
rel_dates = pd.Series(obs_df.index - date0)
# Validate
assert all(rel_dates.dt.seconds==0), "Timestamps are not integer days."
datespan = rel_dates.iloc[-1]-rel_dates.iloc[0]
assert len(rel_dates)-1 == datespan.days, "Some days are missing."
# Set index 
obs_df.index = rel_dates.dt.days
# Crop to start from 0
obs_df = obs_df.loc[0:,:]

# Extract obs
yy = obs_df["deaths"].to_numpy()
# Make non-cumulative
yy = np.diff(yy, prepend=0)
# Obs error matrix
R = 1*eye(1)
Ri = inv(R)
Rm12 = sqrtm(inv(R))
# Obs operator
i_obs = StateVector._fields.index("R_fatl")
def Obs(x):
    return x[...,[i_obs]]

infl = 1.
def analysis(E,y):
    N = len(E)
    N1 = N-1
    for j in [0]:
        Eo = Obs(E)
        xo = mean(Eo,0)
        Y  = Eo - xo
        mu = mean(E,0)
        A  = E-mu
        # Update j-th component of observed ensemble:
        # ------------------------------------------------------
        Y_j    = Rm12[j,:] @ Y.T
        dy_j   = Rm12[j,:] @ (y - xo)
        # Prior var * N1:
        sig2_j = Y_j@Y_j                
        if sig2_j<1e-9: continue
        # Update (below, we drop the locality subscript: _j)
        sig2_u = 1/( 1/sig2_j + 1/N1 )   # KG * N1
        alpha  = (N1/(N1+sig2_j))**(0.5) # Update contraction factor
        dy2    = sig2_u * dy_j/N1        # Mean update
        Y2     = alpha*Y_j               # Anomaly update
        # Without localization:
        Regression = A.T @ Y_j/np.sum(Y_j**2)
        mu        += Regression*dy2
        A         += np.outer(Y2 - Y_j, Regression)

    E = mu + infl*A
    return E



## Integrate
EEf = np.full(tt.shape+E.shape, nan)
EEa = np.full(tt.shape+E.shape, nan)

nVar = len(variables)

for k,t in enumerate(tt[:-1]):
    dt = tt[k+1] - t
    # Variables:
    E[:,:nVar] = rk4(model.dxdt, E[:,:nVar], t, dt)
    # For the params, the model is (as of now) Id.

    # Assimilate
    EEf[k+1] = E
    if k+1<len(yy):
        E = analysis(E,yy[k+1])
        # E.clip(min=1e-9)
    EEa[k+1] = E

    # Write params
    for i,param in enumerate(parameters):
        setattr(model,param,E[:,i+nVar])


EEf[...,:nVar] *= nPop
EEa[...,:nVar] *= nPop

## Plot

# Unpack
state = NamedVars(*EEf.T[:,:60,:])

# Plot
fig, ax = freshfig(1)

coPlot = Lines(ax, state, tt, None, alpha=0.3, lw=1)
# coPlot.add("Exposed")
coPlot.add("Infected")
coPlot.add("Hospitalized")
# coPlot.add("Fatalities",alpha=0.1)
coPlot.add("Recovered")
coPlot.finalize()


t1 = model.t_intervention
t2 = t1 + model.dt_intervention
ax.axvspan(*map(coPlot.t2d, [t1,t2]), alpha=.1, color="b")
ax.text(coPlot.t2d(model.t_intervention), yL[1],"Stricter measures", va="top", ha="right",fontsize="small",rotation=90)
