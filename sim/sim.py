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




def run_simulation(polymer, scatter=None, lines=None, ori_ter=None, show=False, n_steps = 1000,
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
            rnd = random.random()
            if rnd > .5:
                ori_ark = True
            else:
                ori_ark = False
            polymer.coords_tmp = polymer.move_rosenbluth.getOutput(ind1, ind2, cached_counts, ori_ark=ori_ark)
        # print('diff\n ', repr(polymer.coords),  repr(polymer.coords_tmp))

        # rnd = random.random()

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
                # polymer.coords = np.roll(polymer.coords, random.randint(1, polymer.coords.shape[1]), axis=1)
                if show & (step%100 == 0) :
                    # print('show')
                    scatter.x = polymer.coords[0, :];
                    scatter.y = polymer.coords[1, :];
                    scatter.z = polymer.coords[2, :];
                    lines.x = polymer.coords[0, :];
                    lines.y = polymer.coords[1, :];
                    lines.z = polymer.coords[2, :];
                    ori_ter.x = polymer.coords[0, [0, polymer.n //2]]
                    ori_ter.y = polymer.coords[1, [0, polymer.n //2]]
                    ori_ter.z = polymer.coords[2, [0, polymer.n //2]]

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



def make_indexes_and_ark(polymer_):
    """

    :param polymer_:
    :type polymer_:
    :return:
    :rtype:
    """

    ind1 = random.randint(0, polymer_.coords.shape[1] - 1)
    ind2 = ind1
    while ind2 == ind1:
        # print('stuck')
        ind2 = random.randint(0, polymer_.coords.shape[1] - 1)

    # print('before ', a,b, abs(ind2-ind1), ind1, ind2)
    rnd = random.random()
    if rnd > .5:
        ori_ark = True
    else:
        ori_ark = False


    return ind1, ind2, ori_ark


def run_simulation_2(polymer1, polymer2, scatter1=None, lines1=None, scatter2=None, lines2=None, ori_ter=None, show=False, n_steps=1000,
                     use_moves=['move_rosenbluth']):


    """
    running the  simulation
    :return:
    :rtype:
    """
    print('caching  counts...')
    cached_counts = cache_n_conf(polymer1.n, polymer1.cell.A, polymer1.cell.B, polymer1.cell.C)
    print('done')
    print('coords should be equal:', np.array_equal(polymer1.coords, polymer2.coords))
    walks1 = []; walks2 = []
    distances = []
    for step in range(n_steps):
        coords_save1 = polymer1.coords.copy()
        coords_save2 = polymer2.coords.copy()

        if hasattr(polymer1, 'move_rosenbluth') and 'move_rosenbluth' in use_moves:
            polymer1.move_rosenbluth.coordinates = polymer1.coords.copy()

            ind1, ind2, ori_ark = make_indexes_and_ark(polymer1)
           
            polymer1.coords_tmp = polymer1.move_rosenbluth.getOutput(ind1, ind2, cached_counts, ori_ark=ori_ark)

            if step >= n_steps//2:
                # print('dd')
                polymer2.move_rosenbluth.coordinates = polymer2.coords.copy()
                ind1_, ind2_, ori_ark_ = make_indexes_and_ark(polymer2)
                polymer2.coords_tmp = polymer2.move_rosenbluth.getOutput(ind1_, ind2_, cached_counts, ori_ark=ori_ark_)

            # if (ind1 == ind1_) and (ind2 == ind2_):
            #     print(np.array_equal(polymer1.coords_tmp, polymer2.coords_tmp))
            else:
                polymer2.move_rosenbluth.coordinates = polymer1.coords_tmp.copy()
                polymer2.coords_tmp = polymer2.move_rosenbluth.getOutput(0, (polymer2.n)-1, cached_counts, ori_ark=False)


        # print('diff\n ', repr(polymer1.coords),  repr(polymer1.coords_tmp))


        if polymer1.check_borders() &  polymer2.check_borders() & polymer1.check_overlap() :#& polymer2.check_overlap() :
            # calculatin energy for the new configuration
            t1 = np.mean(polymer1.coords_tmp, axis=1)
            dst_new1 = polymer1.cell.f_f.get_distance(t1[1])
            energy_new1 = polymer1.cell.f_f.get_value()(dst_new1)

            t2 = np.mean(polymer2.coords_tmp, axis=1)
            dst_new2 =  polymer2.cell.B  - polymer2.cell.f_f.get_distance(t2[1])
            energy_new2 = polymer2.cell.f_f.get_value()(dst_new2)

            energy_new = energy_new2 + energy_new1

            # calculatin energy for the old configuration
            t1 = np.mean(polymer1.coords, axis=1)
            dst_old1 = polymer1.cell.f_f.get_distance(t1[1])
            energy_old1 = polymer1.cell.f_f.get_value()(dst_old1)

            t2 = np.mean(polymer2.coords, axis=1)
            dst_old2 = polymer2.cell.B  - polymer2.cell.f_f.get_distance(t2[1])
            energy_old2 = polymer2.cell.f_f.get_value()(dst_old2)
            energy_old = energy_old1 + energy_old2

            # print(energy_old, energy_new)
            if (np.exp(-energy_new + energy_old)) >= random.random():
                polymer1.coords = polymer1.coords_tmp.copy()
                polymer2.coords = polymer2.coords_tmp.copy()

                # polymer1.coords = np.roll(polymer1.coords, random.randint(1, polymer1.coords.shape[1]), axis=1)
                if show & (step % 20 == 0):
                    # print(ind1, ind2, ori_ark, ind1_, ind2_, ori_ark_)

                    scatter1.x = polymer1.coords[0, :].copy();
                    scatter1.y = polymer1.coords[1, :].copy();
                    scatter1.z = polymer1.coords[2, :].copy();
                    lines1.x = polymer1.coords[0, :].copy();
                    lines1.y = polymer1.coords[1, :].copy();
                    lines1.z = polymer1.coords[2, :].copy();

                    scatter2.x = polymer2.coords[0, :];
                    scatter2.y = polymer2.coords[1, :];
                    scatter2.z = polymer2.coords[2, :];
                    lines2.x = polymer2.coords[0, :];
                    lines2.y = polymer2.coords[1, :];
                    lines2.z = polymer2.coords[2, :];

                    ori_ter.x = polymer1.coords[0, [0, polymer1.n // 2]]
                    ori_ter.y = polymer1.coords[1, [0, polymer1.n // 2]]
                    ori_ter.z = polymer1.coords[2, [0, polymer1.n // 2]]

                    t1 = np.mean(polymer1.coords, axis=1)
                    t2 = np.mean(polymer2.coords, axis=1)

                    # print("%.1f, %.1f, %.1f"%(t[0], t[1], t[2]))
                    # dst1 = polymer1.cell.f_f.get_distance(t1[1])
                    # dst1 = polymer1.cell.f_f.get_distance(t1[1])

                    walks1.append(t1)
                    walks2.append(t2)

                    # distances.append(dst)
        else:
            # print('outside')
            polymer1.coords = coords_save1
            polymer2.coords = coords_save2


    print("after  moves")
    print(polymer1.coords)
    return walks1, walks2 #distances



def prepare_simulation_2(a,b,c, n):
    """
    prepare simulation
    :return:
    :rtype:
    """

    f_f = ForceField(linear=True, amplitude=17)
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

    rosenbluth1 = Rosenbluth()
    print(rosenbluth1)
    rosenbluth2 = Rosenbluth()


    polymer1 = Polymer(n, cell, kink, crankshaft, pin, rosenbluth1)
    polymer2 = Polymer(n, cell, kink, crankshaft, pin, rosenbluth2)

    # polymer = Polymer(n, cell, rosen)

    print(polymer1)
    print('squashed coordinates are:')
    print(polymer1.get_coords())

    print('unrolled coordinates are:')
    print(polymer1.unroll_chain())

    print(polymer2)
    print('squashed coordinates are:')
    print(polymer2.get_coords())

    print('unrolled coordinates are:')
    print(polymer2.unroll_chain())

    return polymer1, polymer2


if __name__ == "__main__":

    polymer =  prepare_simulation(10,10,10, 50)
    run_simulation(polymer)