##
from corona.utils import *
from corona.maths import *
from corona.model import *
from corona.plotting import *


##

plt.ion()
fig, ax = freshfig(1)


def _iChi2(s,nu):
    return ss.invgamma(a=nu/2, scale=nu/2*s)

def iChi2(mean,var):
    nu = 4 + 2*mean**2/var
    s = mean*(nu-2)/nu
    return _iChi2(s,nu)

# X = ss.lognorm(1.2)
# X = iChi2_from_moments(3,2) 
X = iChi2(0.02, 0.01**2)
xx = X.rvs(10**4)
ax.hist(xx,bins=300)
# print(X.mean(), X.var())

##
def logNorm(mean,var):
    sig2 = log(1 + var/mean**2)
    mu = log(mean/exp(sig2/2))
    return ss.lognorm(s=sqrt(sig2),scale=exp(mu))

X = logNorm(1, .2)
xx = X.rvs(10**4)
ax.hist(xx,bins=300)
print(X.mean(), X.var())

##

##
