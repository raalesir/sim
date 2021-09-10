

## !/usr/bin/env python

"""Tests for `sim` package."""

import pytest
import  numpy as np

from click.testing import CliRunner

from sim.sim import aux
from sim.sim import cli
from sim.sim.consts import N



@pytest.fixture
def response():
    """Sample pytest fixture.

    See more at: http://doc.pytest.org/en/latest/fixture.html
    """
    # import requests
    # return requests.get('https://github.com/audreyr/cookiecutter-pypackage')


def test_content(response):
    """Sample pytest test function with the pytest fixture as an argument."""
    # from bs4 import BeautifulSoup
    # assert 'GitHub' in BeautifulSoup(response.content).title.string


def test_command_line_interface():
    """Test the CLI."""
    runner = CliRunner()
    result = runner.invoke(cli.main)
    assert result.exit_code == 0
    assert 'sim.cli.main' in result.output
    help_result = runner.invoke(cli.main, ['--help'])
    assert help_result.exit_code == 0
    assert '--help  Show this message and exit.' in help_result.output


def test_add():
    assert aux.add(1,2) == 3



def test_check_borders_fail():
    """
    testing ``check_borders``

    :return: True/False
    :rtype: bool
    """

    coords = np.array([[1,1,1], [1,1,2], [1,1,3], [1,1,4], [1,1,5]]).T

    assert aux.check_borders(coords) == False


def test_check_borders_pass():
    """
    testing ``check_borders``

    :return: True/False
    :rtype: bool
    """

    coords = np.array([[1,1,1], [1,1,2], [1,1,3], [1,1,4]]).T

    assert aux.check_borders(coords) == True


def test_check_borders_init_config_pass():
    """
    checking the initial configuration after unrolling

    :return: True/False
    :rtype: bool
    """
    coords = aux.make_circular_chain(N)
    rotation_sequence = aux.generate_x_rotation_seq(N)
    # print(rotation_sequence)
    rotation_matrices = aux.make_rotation_matrices()

    coords = aux.unroll_chain(coords, rotation_sequence, rotation_matrices[:, :, 1])

    print(coords)
    assert aux.check_borders(coords) == True