from corona.utils import *
from scipy.integrate import odeint

# The SIR model differential equations.
def dxdt_param(y, t, beta, gamma):
    S, I, R = y
    N = sum(y)
    dSdt = -beta * S * I / N
    dIdt = beta * S * I / N - gamma * I
    dRdt = gamma * I
    return dSdt, dIdt, dRdt

# Total population, N.
N = 1000

# Initial numbers
I0 = 1           # Infected
R0 = 0           # Recover
S0 = N - I0 - R0 # Susceptible

# Params
beta  = 0.2 # Contact rate
gamma = 0.1 # Mean recovery rate (in people/days).

# A grid of time points (in days)
tt = np.linspace(0, 365, 365+1)

# Fix params
from functools import partial
dxdt = partial(dxdt_param,beta=beta,gamma=gamma)

# Initial conditions vector
x0 = S0, I0, R0

# Integrate the SIR equations over the time grid, t.
ret = odeint(dxdt, x0, tt)
SS, II, RR = ret.T

# Plot the data on three separate curves for S(t), I(t) and R(t)
plt.ion()
fig, ax = freshfig(1)
ax.plot(tt, SS, label='Susceptible')
ax.plot(tt, II, label='Infected')
ax.plot(tt, RR, label='Recovered with immunity')
ax.set_xlabel('Time (days)')
ax.set_ylabel('People')
legend = ax.legend()
