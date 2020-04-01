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

    xx1 = array([6140096.412471, 1262.233441, 818.870653, 5182.640457, 372.438097, 14349.692348, 1268.638699, 663918.896623, 156842.366241, 15887.81097 ])

    assert np.allclose(xx[-1], xx1)
    return xx[-1]

