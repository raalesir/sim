"""
    testing  routines  in aux
"""

import numpy as np

try:
    from sim.sim import aux
except ModuleNotFoundError:
    try:
        from sim import aux
    except:
        from .sim  import  aux



def test_cache_n_conf():
    res1 = np.array(
        [[[[0., 1.],
           [1., 0.]],
          [[1., 0.],
           [0., 0.]]]]
    )
    res2 = aux.cache_n_conf(1, 2, 2, 2)
    assert  np.array_equal(res1, res2)==True



def test_cache_n_conf2():
    res1 = np.array(
        [[[6, 0],
          [0, 2]],
         [[0, 2],
          [2, 0]]]
    )

    res2 = aux.cache_n_conf(10, 2, 2, 2)[1, :, :, :]

    assert np.array_equal(res1, res2) ==True


def test_cache_n_conf1():
    res2 = aux.cache_n_conf(10, 5, 5, 5)
    assert res2.shape == (10, 5, 5, 5)




def test_n_conf_pass1():
    """ ring polymer,  analytical result   is known"""

    assert aux.n_conf(6, 0, 0, 0) == 1860.0



def test_n_conf_pass2():
    """
    stretched conformation, only one conf possible
    """
    assert aux.n_conf(3, 3 ,0,0) == 1



def test_n_conf_pass3():
    """
    not possible conf, zero conformations
    """
    assert aux.n_conf(3, 2 ,0,0) == 0.0


def test_n_conf_pass4():
    """
    checked on paper
    """
    assert aux.n_conf(4, 2 ,0,0) == 28.0


def test_n_conf_pass5():
    """
    checked on paper
    """
    assert aux.n_conf(4, 0 ,0,2) == 28.0