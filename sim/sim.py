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
    from sim.consts import N
    from sim.cell import CubicCell
    # from sim.aux import *
    from sim.moves import Kink, CrankShaft, Pin
    from sim.polymer import Polymer

except ModuleNotFoundError:
    from cell import CubicCell

    from consts import  N
    # from aux import *
    from moves import Kink,CrankShaft, Pin
    from polymer import Polymer




def run_simulation():
    """
    running the  simulation
    :return:
    :rtype:
    """
    cell = CubicCell(6, 10, 6)
    print(cell)
    print(cell.get_center())

    kink = Kink()
    print(kink)
    crankshaft = CrankShaft()
    print(crankshaft)

    pin = Pin()
    print(pin)

    polymer = Polymer(N, cell, kink, crankshaft, pin)
    print(polymer)
    print('squashed coordinates are:')
    print(polymer.get_coords())

    print('unrolled coordinates are:')
    print(polymer.unroll_chain())


    for step in range(1000):
        coords_save = polymer.coords.copy()
        rnd = random.random()
        if rnd < 0.4:
            crankshaft.coordinates = polymer.coords
            polymer.coords_tmp = polymer.move_crankshaft.getOutput(1)
        elif (rnd > 0.4) & (rnd < 0.7):
            kink.coordinates = polymer.coords  #
            polymer.coords_tmp = polymer.move_kink.getOutput(1)
        else:
            #         print('pin')
            pin.coordinates = polymer.coords
            polymer.coords_tmp = polymer.move_pin.getOutput(1)

        if polymer.check_borders() & polymer.check_overlap():
            polymer.coords = polymer.coords_tmp.copy()

        else:
            polymer.coords = coords_save

    print("after  moves")
    print(polymer.coords)



if __name__ == "__main__":

    run_simulation()