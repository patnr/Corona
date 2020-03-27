import numpy as np
from scipy.integrate import odeint
import matplotlib.pyplot as plt
from corona.utils import freshfig

def SEIR():
    Rt = 2
    Tinc = 5
    Tinf = 3
    N = 5000001
    def dxdt(x, t, *args):
        S, E, I, R = x
        Rt, Tinf, Tinc , N = args
       
        dSdt = -Rt*I*S/(Tinf*N)
        dEdt = Rt*I*S/(Tinf*N)-E/Tinc
        dIdt = E/Tinc - I/Tinf
        dRdt = I/Tinf
        return dSdt, dEdt, dIdt, dRdt

    t = np.linspace(0, 365, 1000) # start, stop, numSamples
    S0 = 5000000
    E0 = 0
    I0 = 1
    R0 = 0
    x0 = S0, E0, I0, R0
    N = S0 + E0 + I0 + R0
    sol = odeint(dxdt, x0, t, args=( Rt, Tinf, Tinc, N ))

    ax.plot(t, sol[:, 0], 'r', label='S(t)')
    ax.plot(t, sol[:, 1], 'g', label='E(t)')
    ax.plot(t, sol[:, 2], 'b', label='I(t)')
    ax.plot(t, sol[:, 3], 'k', label='R(t)')
    ax.plot(t, sol[:, 0] + sol[:,1] + sol[:, 2] + sol[:, 3], 'm', label='population')
    ax.legend(loc='best')
    ax.grid()
    ax.set_xlabel('t')

plt.ion()
fig, ax = freshfig(2)
SEIR()

