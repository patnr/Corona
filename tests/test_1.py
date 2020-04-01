# -*- coding: utf-8 -*-

import pytest

from corona.utils import *
from corona.model import *

__author__ = "patricknraanes"
__copyright__ = "patricknraanes"
__license__ = "mit"

## Params
def test_SEIR2():
    model = SEIR2()

    # Population size
    nPop = 7*10**6

    # Time -- unit: days
    t_end = 200
    dt = 0.1
    tt = linspace(0, t_end, int(t_end/dt)+1)

    ## Integrate
    x0 = model.init_state(Infected=1/nPop)
    xx = integrate(model.dxdt, x0, tt)

    # Multiply by population
    xx = nPop * xx

    xx1 = array([6148274.868637, 1257.648524, 815.541392, 5147.231294, 370.609155, 14203.185958, 1255.869424, 657581.283035, 155356.593176, 15737.169405])

    assert np.allclose(xx[-1], xx1)
    # return xx[-1]

