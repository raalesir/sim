"""
Main module
===========
"""

# import ipyvolume as ipv
import numpy as np

# from time  import  sleep

# import matplotlib.pyplot as plt

import random


try:
    # from sim.consts import N
    from sim.cell import CubicCell, ForceCubicCell
    from sim.aux import  cache_n_conf
    from sim.moves import Kink, CrankShaft, Pin, Rosenbluth
    from sim.polymer import Polymer
    from sim.force_field import  ForceField

except ModuleNotFoundError:
    try:
        from cell import CubicCell, ForceCubicCell
        # from consts import  N
        from aux import cache_n_conf
        from moves import Kink,CrankShaft, Pin, Rosenbluth
        from polymer import Polymer
        from force_field import  ForceField
    except:
        from .cell import  CubicCell, ForceCubicCell
        from .moves import  Kink, CrankShaft,Pin, Rosenbluth
        from .polymer import Polymer
        from .aux import  cache_n_conf
        from .force_field import ForceField



def prepare_simulation(a,b,c, n):
    """
    prepare simulation
    :return:
    :rtype:
    """

    f_f = ForceField(linear=True, amplitude=20)
    f_f.origin = b
    print(f_f)


    # cell = CubicCell(a, b, c)
    cell = ForceCubicCell(a,b,c, f_f)
    print(cell)
    print(cell.get_center())
    print("the force origin is at: %s" %cell.f_f.origin)


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




def run_simulation(polymer, scatter=None, lines=None, show=False, n_steps = 1000,
                   use_moves=['move_rosenbluth']):
    """
    running the  simulation
    :return:
    :rtype:
    """
    print('caching  counts...')
    cached_counts = cache_n_conf(polymer.n, polymer.cell.A, polymer.cell.B, polymer.cell.C)
    print('done')

    walks=[]
    distances=[]
    for step in range(n_steps):
        coords_save = polymer.coords.copy()
        if hasattr(polymer, 'move_rosenbluth') and 'move_rosenbluth' in use_moves:
            polymer.move_rosenbluth.coordinates = polymer.coords.copy()

            ind1 = random.randint(0, polymer.coords.shape[1]-1)
            ind2 = ind1
            while ind2 == ind1:
                # print('stuck')
                ind2 = random.randint(0, polymer.coords.shape[1]-1)

            # print('before ', a,b, abs(ind2-ind1), ind1, ind2)
            polymer.coords_tmp = polymer.move_rosenbluth.getOutput(ind1, ind2, cached_counts)
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
            # calculatin energy for the new configuration
            t = np.mean(polymer.coords_tmp, axis=1)
            dst_new1 = polymer.cell.f_f.get_distance(t[1])
            energy_new1 = polymer.cell.f_f.get_value()(dst_new1)
            dst_new2 = polymer.cell.B - dst_new1
            energy_new2 = polymer.cell.f_f.get_value()(dst_new2)
            energy_new = energy_new2 + energy_new1

            # calculatin energy for the old configuration
            t = np.mean(polymer.coords, axis=1)
            dst_old1 = polymer.cell.f_f.get_distance(t[1])
            energy_old1 = polymer.cell.f_f.get_value()(dst_old1)
            dst_old2 = polymer.cell.B - dst_old1
            # print(dst_old2, t[1])
            energy_old2 = polymer.cell.f_f.get_value()(dst_old2)
            energy_old = energy_old1 + energy_old2

            # print(energy_old, energy_new)
            if (energy_new-energy_old) < rnd:
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
                    dst  = polymer.cell.f_f.get_distance(t[1])
                    walks.append(t)
                    distances.append(dst)
        else:
            # print('outside')
            polymer.coords = coords_save

    print("after  moves")
    print(polymer.coords)
    return walks, distances



if __name__ == "__main__":

    polymer =  prepare_simulation(10,10,10, 50)
    run_simulation(polymer)