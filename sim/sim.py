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
    from sim.aux import  cache_n_conf
    from sim.moves import Kink, CrankShaft, Pin, Rosenbluth
    from sim.polymer import Polymer

except ModuleNotFoundError:
    try:
        from cell import CubicCell
        # from consts import  N
        from aux import cache_n_conf
        from moves import Kink,CrankShaft, Pin, Rosenbluth
        from polymer import Polymer
    except:
        from .cell import  CubicCell
        from .moves import  Kink, CrankShaft,Pin, Rosenbluth
        from .polymer import Polymer
        from .aux import  cache_n_conf



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

    rosenbluth = Rosenbluth()
    print(rosenbluth)


    polymer = Polymer(n, cell, kink, crankshaft, pin, rosenbluth)
    # polymer = Polymer(n, cell, rosen)

    print(polymer)
    print('squashed coordinates are:')
    print(polymer.get_coords())

    print('unrolled coordinates are:')
    print(polymer.unroll_chain())

    return polymer




def run_simulation(polymer, scatter=None, lines=None, show=False, n_steps = 1000, use_moves=['move_rosenbluth']):
    """
    running the  simulation
    :return:
    :rtype:
    """
    print('caching  counts...')
    cached_counts = cache_n_conf(polymer.n, polymer.cell.A, polymer.cell.B, polymer.cell.C)
    print('done')

    walks=[]
    for step in range(n_steps):
        coords_save = polymer.coords.copy()
        if hasattr(polymer, 'move_rosenbluth') and 'move_rosenbluth' in use_moves:
            polymer.move_rosenbluth.coordinates = polymer.coords.copy()

            ind1 = random.randint(0, polymer.coords.shape[1]-1)
            ind2 = ind1
            while ind2 == ind1:
                # print('stuck')
                ind2 = random.randint(0, polymer.coords.shape[1]-1)

            a = int(polymer.coords[0, ind1]), int(polymer.coords[1, ind1]), int(polymer.coords[2, ind1])
            b = int(polymer.coords[0, ind2]), int(polymer.coords[1, ind2]), int(polymer.coords[2, ind2])
            # print('before ', a,b, abs(ind2-ind1), ind1, ind2)
            polymer.coords_tmp = polymer.move_rosenbluth.getOutput(a,b, ind1, ind2, cached_counts)
        # print('diff\n ', repr(polymer.coords),  repr(polymer.coords_tmp))

        rnd = random.random()

        # if rnd < 0.4:
        #     # if hasattr(polymer, 'move_crankshaft'):
        #     polymer.move_crankshaft.coordinates = polymer.coords
        #         # crankshaft.coordinates = polymer.coords
        #     polymer.coords_tmp = polymer.move_crankshaft.getOutput(n_steps=1)
        # elif (rnd > 0.4) & (rnd < 0.7):
        #     polymer.move_kink.coordinates =  polymer.coords
        #     # kink.coordinates = polymer.coords  #
        #     polymer.coords_tmp = polymer.move_kink.getOutput(n_steps=1)
        # else:
        #     # print('pin')
        #     polymer.move_pin.coordinates  = polymer.coords
        #     # pin.coordinates = polymer.coords
        #     polymer.coords_tmp = polymer.move_pin.getOutput(n_steps =1)

        if polymer.check_borders() & polymer.check_overlap():
            polymer.coords = polymer.coords_tmp.copy()
            polymer.coords = np.roll(polymer.coords, random.randint(1, polymer.coords.shape[1]), axis=1)
            if show & (step%100 == 0) :
                # print('show')
                scatter.x = polymer.coords[0, :];
                scatter.y = polymer.coords[1, :];
                scatter.z = polymer.coords[2, :];
                lines.x = polymer.coords[0, :];
                lines.y = polymer.coords[1, :];
                lines.z = polymer.coords[2, :];
                t =  np.mean(polymer.coords, axis=1)
                # print("%.1f, %.1f, %.1f"%(t[0], t[1], t[2]))
                walks.append(t)
        else:
            # print('outside')
            polymer.coords = coords_save

    print("after  moves")
    print(polymer.coords)
    return walks



if __name__ == "__main__":

    polymer =  prepare_simulation(10,10,10, 50)
    run_simulation(polymer)