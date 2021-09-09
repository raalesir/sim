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

