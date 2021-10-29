"""
    Testing moves
"""

import pytest
import numpy as np


try:
    from sim.sim import moves,consts, cell, polymer
except ModuleNotFoundError:
    try:
        from sim import moves, cell, consts, polymer
    except:
        from .sim  import  moves, cell, consts, polymer


def test_kink_move_pass():
    """
    testing kink move

    :return: True/False
    :rtype: bool
    """

    kink = moves.Kink()

    kink.coordinates = np.array([[0,0,1], [1,0,0], [0,0,0]])

    res = kink.getOutput(1)

    assert np.array_equal(res, np.array([[0,1,1], [1,1,0], [0,0,0]]))


def test_kink_move_fail():
    """
    testing kink move

    :return: True/False
    :rtype: bool
    """

    kink = moves.Kink()

    kink.coordinates = np.array([[0,0,1], [1,0,0], [0,0,0]])

    res = kink.getOutput(1)

    assert not  np.array_equal(res, np.array([[0,0,1], [1,0,0], [0,0,0]]))


def test_crank_move_pass():
    """
    checking the crankshaft move

    :return: True/False
    :rtype: bool
    """
    crankshaft = moves.CrankShaft()

    crankshaft.coordinates  = np.array([[0,0,0,0,0,0],
                                       [0,1,1,2,2,3],
                                       [0,0,1,1,0,0]])
    res = crankshaft.crankshaft_move(1,4, consts.rot[:,:,3])

    res1 = np.array([[0,0,-1,-1,0,0],
                    [0,1,1,2,2,3],
                    [0,0,0,0,0,0]])
    assert np.array_equal(res, res1)



def test_crank_move_pass1():
    """
    checking the crankshaft move

    :return: True/False
    :rtype: bool
    """
    crankshaft = moves.CrankShaft()

    crankshaft.coordinates = np.array([[0, 0, 0, 0, 0, 0],
                                       [0, 1, 1, 2, 2, 3],
                                       [0, 0, 1, 1, 0, 0]])
    res = crankshaft.crankshaft_move(0, 5, consts.rot[:, :, 3])

    res1 = np.array([[0, 0, -1, -1, 0, 0],
                     [0, 1, 1, 2, 2, 3],
                     [0, 0, 0, 0, 0, 0]])
    assert np.array_equal(res, res1)



def test_crank_move_pass1():
    """
    checking the crankshaft move

    :return: True/False
    :rtype: bool
    """
    crankshaft = moves.CrankShaft()

    crankshaft.coordinates = np.array([[0, 0, 0, 0, 0, 0],
                                       [0, 1, 1, 2, 2, 3],
                                       [0, 0, 1, 1, 0, 0]])
    res = crankshaft.crankshaft_move(0, 5, consts.rot[:, :, 3])

    res1 = np.array([[0, 0, -1, -1, 0, 0],
                     [0, 1, 1, 2, 2, 3],
                     [0, 0, 0, 0, 0, 0]])
    assert np.array_equal(res, res1)



def test_poolymer_check_borders_fail():
    """
    testing ``check_borders``

    :return: True/False
    :rtype: bool
    """

    cell1 =  cell.CubicCell(2,2,2)
    polymer1 = polymer.Polymer(5, cell1)
    polymer1.coords_tmp = np.array([[0, 1,1], [1,1,1], [1,1,2], [1,1,3], [1,1,4]]).T
    res = polymer1.check_borders()

    assert res == False



def test_poolymer_check_borders_pass():
    """
    testing ``check_borders``

    :return: True/False
    :rtype: bool
    """

    cell1 =  cell.CubicCell(5,5,5)
    polymer1 = polymer.Polymer(5, cell1)
    polymer1.coords_tmp = np.array([[1, 1,1], [1,1,2], [1,1,3], [1,2,3], [1,3,3]]).T
    res = polymer1.check_borders()

    assert res == True

