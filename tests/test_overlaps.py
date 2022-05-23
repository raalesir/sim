"""
    testing for the overlap module
"""

# import pytest

try:
    from sim.sim import overlaps
except ModuleNotFoundError:
    from  sim import   overlaps


def test_n_conf():
    """
    testing number of conformations for a given `N`
    :return: True/False
    :rtype: bool
    """

    overlaps_10 = overlaps.Overlap(N=10)

    assert  overlaps_10.n_conformations == 1172556


def test_overlap_distribution():
    """
    checking the distribution for overlaps
    :return: True/False
    :rtype: bool
    """

    overlaps_4 = overlaps.Overlap(N=4)

    assert  overlaps_4.get_overlaps_histogram() == {0: 24, 1: 60, 2: 6}



def test_n_conf_all_confs():
    """
    the number of confs for phantom polymer according to the  formula should be equal to the length of the
    list with all  confs  constructed

    :return: True/False
    :rtype: bool
    """

    overlaps_4 = overlaps.Overlap(N=4)
    overlaps_4.calculate_all_conformations()
    assert len(overlaps_4.keep_result.all) == overlaps_4.n_conformations

