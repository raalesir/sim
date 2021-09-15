"""
    Moves class
    ===========
"""

import numpy as np
import random

try:
    from sim.consts  import rot
except ModuleNotFoundError:
    try:
        from consts import rot
    except ModuleNotFoundError:
        from  .consts import rot


class Move:
    """
    basic move class
    """

    def __init__(self):
        """

        :param coordinates: coordinates to change
        :type coordinates: numpy array, int
        """
        self.coordinates = None
        self.output = None
        self.length = None


    def getOutput(self, n_steps=100):
        """
        performing a move

        :return: 3D coordinates
        :rtype: numpy array, int
        """

        self.length = self.coordinates.shape[1]

        self.output = self.make_move(n_steps)
        return  self.output



class  CrankShaft(Move):
    """
    class for crankshaft move
    """

    def __init__(self):
        super(Move, self).__init__()

    def __str__(self):
        return "crankshaft"


    def crankshaft_move(self, A, B, rotation):
        """
            performs a crankshaft move for given part of the chain around given axis on given angle.

        :param coords_:    (3,N) coordinates
        :type coords_:      numpy array
        :param A:         starting position to  crank  about
        :type A:             int
        :param B:         ending position to crank about
        :type B:            int
        :param rotation:    (3,3) rotation matrix
        :type rotation:      Numpu array
        :return:            (3,N) coordinates after crankshaft
        :rtype:              Numpy array
        """

        tail = self.coordinates[:, A + 1:B]
        tail0 = self.coordinates[:, A:A + 1]
        #     print(tail, tail0)
        tail = tail - tail0
        #     print(tail)

        #     print(rotation.shape)
        for i in range(tail.shape[1]):
            self.coordinates[:, A + i + 1] = np.matmul(rotation, tail[:, i]) + tail0[:, 0]

        return self.coordinates



    def prepare_crank(self):
        """
           generates input data for crankshaft move --  crank points and the index in the rotation list


        :return:    tuple of 3 integers, ``(crank, crank1, rotation)``, where ``crank`` is the starting and ``crank1`` the ending \
        indexes of the coordinates array for the crankshaft move: ``[crank+1:crank1-1]`` are subject to crankshaft; ``rotation`` -- \
        is the index of the list with rotation matrices
        :rtype:  tuple
        """

        potentials = []  # if this list is empty, it means that it the  chosen bead `crank` has to buddy
        # to form the rotation axis
        while not potentials:
            # selecting random bead
            crank = np.random.randint(self.length)
            #     print(crank)

            # generating potential positions for  the second bead
            potentials_x = np.where((self.coordinates[1] == self.coordinates[1, crank]) & (self.coordinates[2] == self.coordinates[2, crank]))[0]
            potentials_y = np.where((self.coordinates[0] == self.coordinates[0, crank]) & (self.coordinates[2] == self.coordinates[2, crank]))[0]
            potentials_z = np.where((self.coordinates[0] == self.coordinates[0, crank]) & (self.coordinates[1] == self.coordinates[1, crank]))[0]

            potentials = [el for el in [(potentials_x, [0, 1]),
                                        (potentials_y, [2, 3]),
                                        (potentials_z, [4, 5])] if len(el[0]) > 0
                          ]

            # selecting the axis with rotatios
            axis_rotation = random.choice(potentials)
            crank1 = random.choice(axis_rotation[0])
            rotation = random.choice(axis_rotation[1])  # index in the list of rotation matrices

        return crank, crank1, rotation



    def make_move(self, n_steps):
        """

        :param c: (3,N) coords for ring polymer
        :type c: numpy array, int
        :param n_steps: number of crankshafts
        :type n_steps: int
        :param rot: rotation matrices
        :type rot: numpy (3,3,6)  array, float
        :return: (3,N) coords after performing ``n_steps`` random crankshafts
        :rtype: numpy array, int
        """

        # rot = self._make_rotation_matrices()

        for i in range(n_steps):
            crank, crank1, rotation = self.prepare_crank()
            self.coordinates = self.crankshaft_move(crank, crank1, rot[:, :, rotation])

        return self.coordinates



class Kink(Move):
    """
    class for kink move
    """

    def __init__(self):
        super(Move, self).__init__()


    def __str__(self):
        return "kink"


    def make_move(self, n_steps=100):
        """
            runs the ``n_steps`` random kinks

        :param n_steps: number of kinks
        :type n_steps: int
        :return: (3,N) coords after performing ``n_steps`` random kinks
        :rtype:   numpy array, int
        """



        for i in range(n_steps):
            kink, kink_minus, kink_plus = self.prepare_kink()

            self.coordinates = self.kink_move(kink, kink_minus, kink_plus)

        return self.coordinates



    def kink_move(self, kink, kink_minus, kink_plus):
        """
                performs a kink move


        :param kink: the bead subjected to kink move
        :type kink: int
        :param kink_minus: preceding bead
        :type kink_minus: int
        :param kink_plus: succeeding bead
        :type kink_plus: int
        :return: (3,N) coordinates after the kink move
        :rtype: Numpy array
        """

        if self.coordinates[0, kink_minus] == self.coordinates[0, kink] == self.coordinates[0, kink_plus]:
            self.coordinates[1, kink] = self.coordinates[1, kink_plus]
            self.coordinates[2, kink] = self.coordinates[2, kink_minus]
        # print('x')

        elif self.coordinates[1, kink_minus] == self.coordinates[1, kink] == self.coordinates[1, kink_plus]:
            self.coordinates[0, kink] = self.coordinates[0, kink_plus]
            self.coordinates[2, kink] = self.coordinates[2, kink_minus]
        # print('y')


        else:
            self.coordinates[0, kink] = self.coordinates[0, kink_plus]
            self.coordinates[1, kink] = self.coordinates[1, kink_minus]
        # print('z')

        return self.coordinates



    def prepare_kink(self):
        """
            prepares input data for the kink  move

        :return: 3-tuple as (preceding, actual, succeeding)
        :rtype: tuple of int
        """

        orthogonal = False

        while not orthogonal:

            kink = np.random.randint(self.length)

            kink_minus = kink - 1
            kink_plus = kink + 1

            if kink == 0:
                kink_minus = self.length - 1
            elif kink == self.length - 1:
                kink_plus = 0

                #     print(kink, length, coords[:, kink])

            v1 = self.coordinates[:, kink] - self.coordinates[:, kink_minus]
            v2 = self.coordinates[:, kink_plus] - self.coordinates[:, kink]

            orthogonal = np.dot(v1, v2) == 0

        return kink, kink_minus, kink_plus

