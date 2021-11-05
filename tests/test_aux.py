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



def test_get_n_beads1():

    assert aux.get_n_beads(0, 0, 0) == 1


def test_get_n_beads2():
    assert aux.get_n_beads(10, 0, 0) == 1

def test_get_n_beads3():
    assert aux.get_n_beads(10, 0, 5) == 6

def test_get_n_beads5():
    assert aux.get_n_beads(10, 1, 5) == 5

def test_get_n_beads6():
    assert aux.get_n_beads(10, 5, 1) == 5


def test_get_n_beads4():
    "9,0,1,2"
    assert aux.get_n_beads(10, 2, 9) == 4

def test_get_n_beads7():
    "9,0,1,2"
    assert aux.get_n_beads(10, 9, 2) == 4


def test_get_n_beads8():
    "9,0,1,2"
    assert aux.get_n_beads(10, 12, 9) == 0

def test_get_n_beads9():
    "9,0,1,2"
    assert aux.get_n_beads(10, -1, 9) == 0



def test_get_sequence_of_coords():

    assert aux.get_sequence_of_coords(0,0,0) == [0]


def test_get_sequence_of_coords1():
    assert aux.get_sequence_of_coords(10, -1, 1)  == None


def test_get_sequence_of_coord2():
    assert aux.get_sequence_of_coords(10, 0, 4, ori_ark=False) == [0,1,2,3,4]


def test_get_sequence_of_coord3():
    assert aux.get_sequence_of_coords(10, 0, 5, ori_ark=True) == [5, 6, 7, 8, 9, 0]



def test_get_sequence_of_coord4():
    assert aux.get_sequence_of_coords(10, 1, 4, ori_ark=True) == [4, 5, 6, 7, 8, 9, 0, 1]


def test_get_sequence_of_coord5():
    assert aux.get_sequence_of_coords(10, 4, 0, ori_ark=True) == [4, 5, 6, 7, 8, 9, 0]



def test_get_ind1ind2_1():

    assert aux.get_ind1ind2(10,1000, 100) == (99,0)


def test_get_ind1ind2_2():

    assert aux.get_ind1ind2(0,1000, 100) == (99,0)


def test_get_ind1ind2_3():
    assert aux.get_ind1ind2(500, 1000, 100) == (74, 25)


def test_get_ind1ind2_4():
    assert aux.get_ind1ind2(1000, 1000, 100) == (49, 50)
