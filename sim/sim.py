"""
Main module
===========
"""

# import ipyvolume as ipv
import numpy as np

from time  import  sleep

import matplotlib.pyplot as plt

import random


try:
    # from sim.consts import N
    from sim.cell import CubicCell
    # from sim.aux import *
    from sim.moves import Kink, CrankShaft, Pin
    from sim.polymer import Polymer

except ModuleNotFoundError:
    try:
        from cell import CubicCell
        # from consts import  N
        # from aux import *
        from moves import Kink,CrankShaft, Pin
        from polymer import Polymer
    except:
        from .cell import  CubicCell
        from .moves import  Kink, CrankShaft,Pin
        from .polymer import Polymer


def prepare_simulation(a,b,c, n):
    """
    prepare simulation
    :return:
    :rtype:
    """
    cell = CubicCell(a, b, c)
    print(cell)
    print(cell.get_center())

    kink = Kink()
    print(kink)
    crankshaft = CrankShaft()
    print(crankshaft)
    pin = Pin()
    print(pin)

    polymer = Polymer(n, cell, kink, crankshaft, pin)
    print(polymer)
    print('squashed coordinates are:')
    print(polymer.get_coords())

    print('unrolled coordinates are:')
    print(polymer.unroll_chain())

    return polymer



def run_simulation(polymer, scatter=None, lines=None, show=False, n_steps = 1000) :
    """
    running the  simulation
    :return:
    :rtype:
    """


    for step in range(n_steps):
        coords_save = polymer.coords.copy()
        rnd = random.random()
        if rnd < 0.4:
            polymer.move_crankshaft.coordinates = polymer.coords
            # crankshaft.coordinates = polymer.coords
            polymer.coords_tmp = polymer.move_crankshaft.getOutput(1)
        elif (rnd > 0.4) & (rnd < 0.7):
            polymer.move_kink.coordinates =  polymer.coords
            # kink.coordinates = polymer.coords  #
            polymer.coords_tmp = polymer.move_kink.getOutput(1)
        else:
            # print('pin')
            polymer.move_pin.coordinates  = polymer.coords
            # pin.coordinates = polymer.coords
            polymer.coords_tmp = polymer.move_pin.getOutput(1)

        if polymer.check_borders() & polymer.check_overlap():
            polymer.coords = polymer.coords_tmp.copy()
            if show:
                scatter.x = polymer.coords[0, :];
                scatter.y = polymer.coords[1, :];
                scatter.z = polymer.coords[2, :];
                lines.x = polymer.coords[0, :];
                lines.y = polymer.coords[1, :];
                lines.z = polymer.coords[2, :];

        else:
            polymer.coords = coords_save

    print("after  moves")
    print(polymer.coords)



if __name__ == "__main__":

    polymer =  prepare_simulation(10,10,10, 50)
    run_simulation(polymer)