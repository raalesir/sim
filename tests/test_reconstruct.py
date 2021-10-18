
"""
    testing for reconstruction module
"""

import numpy as np
try:
    from sim.sim import aux
except ModuleNotFoundError:
    from  sim import   aux


def test_loss_function_pass():

    t1 = np.diag((3,3,3))
    r = aux.loss_function( (1, 0, 0, 0, 1, 0, 0, 0, 1,0,0,0), t1, t1)

    assert r  == 0.0


def test_loss_function_pass1():

    t1 = np.diag((3, 3, 3))
    t2 = np.diag((1, 1, 1))

    r = aux.loss_function((1, 0, 0, 0, 1, 0, 0, 0, 1,0,0,0), t1, t2)
    assert r == 2.0


def test_loss_function_fail():

    t1 = np.diag((3, 3, 3))
    t2 = np.diag((1, 1, 1))

    r = aux.loss_function((1, 0, 0, 0, 1, 0, 0, 0, 1,0,0,0), t1, t2)
    assert r != 1.0


