"""
    Testing moves
"""

import pytest
import numpy as np


try:
    from sim import moves
except ModuleNotFoundError:
    from sim import moves


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


