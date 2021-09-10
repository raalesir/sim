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
    from sim.aux import *

except ModuleNotFoundError:
    from consts import  N
    from aux import *

if __name__ == "__main__":

    print(N)

    coords = make_circular_chain(N)
    rotation_sequence = generate_x_rotation_seq(N)
    print(rotation_sequence)

    rotation_matrices = make_rotation_matrices()

    coords = unroll_chain(coords, rotation_sequence, rotation_matrices[:, :, 1])

    print(coords)

    print('kinking')
    coords = run_kinks(coords, 10000)
    print(coords)

    # print("running crankshafts")
    # coords  = run_crankshafts(coords, 10, rotation_matrices)
    # print(coords)
    print("checking borders")
    print(check_borders(coords))

    print("running crankshafts")

    for  step in range(10000):
        coords_save = coords.copy()
        if random.random() < 0.2:
            coords_tmp = run_crankshafts(coords, 1, rotation_matrices)
        else:
            coords_tmp = run_kinks(coords, 1)

        if check_borders(coords_tmp):
            coords = coords_tmp
            # print('post check', check_borders(coords))

        else:
            coords = coords_save

    print("checking borders")
    print(check_borders(coords))


