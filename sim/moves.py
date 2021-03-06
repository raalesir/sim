"""
    Moves class
    ===========
"""

import numpy as np
import random
import sys

try:
    from sim.consts  import rot
    from  sim.aux import n_conf, get_n_beads, get_sequence_of_coords
    # from sim import aux

except ModuleNotFoundError:
    try:
        from consts import rot
        # import  aux
        from  aux import n_conf, get_n_beads, get_sequence_of_coords
    except ModuleNotFoundError:
        from  .consts import rot
        from .aux import n_conf, get_n_beads, get_sequence_of_coords
        # import  .aux



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


    # def getOutput(self, n_steps=100,):
    def getOutput(self, *args, **kwargs):

        """
        performing a move

        :return: 3D coordinates
        :rtype: numpy array, int
        """
        if 'n_steps' in kwargs:
            n_steps = kwargs['n_steps']

        else:
            n_steps=1

        self.length = self.coordinates.shape[1]
        self.output = self.make_move(n_steps)
        return  self.output



class Rosenbluth(Move):
    """
    class  for regrowing chain between two points
    """
    use_cache = True

    def __init__(self):
        super(Move, self).__init__()


    def __str__(self):
        return "rosenbluth"


    @staticmethod
    def get_neighbours(n, dx, dy, dz):
        """
        returns all possible  next steps
        """
        r = []
        r.append([n - 1, dx + 1, dy, dz])
        r.append([n - 1, dx - 1, dy, dz])

        r.append([n - 1, dx, dy + 1, dz])
        r.append([n - 1, dx, dy - 1, dz])

        r.append([n - 1, dx, dy, dz + 1])
        r.append([n - 1, dx, dy, dz - 1])

        return r


    @staticmethod
    def make_limits(a_, b_):

        k1, l1, m1 = a_
        k2, l2, m2 = b_

        return k2 - k1, l2 - l1, m2 - m1

    @staticmethod
    def build_chain(N, dx, dy, dz, coords, counts):
        """
        build chain between 2 ends
        """

        if N < 2:
            #         print('base case')
            #         print("0,0,0,0")
            return [[0, 0, 0]]
            #         return

        else:
            rr = []
            neighbours = Rosenbluth.get_neighbours(N, dx, dy, dz)
            for neighbour in neighbours:
                if Rosenbluth.use_cache:
                    try:
                        t = counts[neighbour[0] - 1, abs(neighbour[1]), abs(neighbour[2]), abs(neighbour[3])]
                    except:
                        t = n_conf(*neighbour)
                        print('miss. out of box? ', neighbour)
                    rr.append(t)
                else:
                    rr.append(n_conf(*neighbour))
            try:
                l = [item / sum(rr) for item in rr]
            except:
                l = [0 * len(rr)]
            r = np.cumsum(l)
            try:
                selected = neighbours[np.where(r > np.random.random())[0][0]]
                c = [selected[1:]] + Rosenbluth.build_chain(*selected, selected[1:], counts)
            except:
                # no way to build chain
                print('no way to  build chain')
                selected = [1, 0, 0, 0]
                c = [selected[1:]] + Rosenbluth.build_chain(*selected, selected[1:], counts)

            return c


    def getOutput(self, i1, i2, cached_counts):
        """
        regrows the chain betweeen the i1 and i2

        :param i1: start index  to regrow
        :type i1: int
        :param i2: last index to regrow
        :type i2: int
        :param cached_counts: Array with  cached number  of trajectories for given `N, k,l,m`
        :type cached_counts:  4D numpy array
        :return: regrown coordinates
        :rtype: (3, N) Numpy array of integers
        """

        self.length = self.coordinates.shape[1]
        n_beads = get_n_beads(self.length, i1, i2) - 1
        # n_beads = abs(i1 -i2)
        # print(a,b,i1,i2, 'n_beads=', n_beads)
        # print(self.coordinates)

        a = int(self.coordinates[0, i1]), int(self.coordinates[1, i1]), int(self.coordinates[2, i1])
        b = int(self.coordinates[0, i2]), int(self.coordinates[1, i2]), int(self.coordinates[2, i2])


        # def reverse(A):
        #     """
        #     reverses the trajectory with mirroring
        #     """
        #     diffs = np.diff(A, 1, axis=0)
        #     rev = [A[-1]]
        #
        #     for d in diffs:
        #         rev.append(rev[-1] - d)
        #
        #     return np.array(rev)[::-1]


        dk, dl, dm = Rosenbluth.make_limits(a, b)

        regrown_coords = np.array([[abs(dk), abs(dl), abs(dm)]] + \
                               Rosenbluth.build_chain(n_beads, abs(dk), abs(dl), abs(dm), [], cached_counts)
                               )

        # print(regrown_coords.shape, n_beads, abs(i1-i2))
        if dk < 0:
            regrown_coords[:, 0] = -  regrown_coords[:, 0]
        if dl < 0:
            regrown_coords[:, 1] = -  regrown_coords[:, 1]
        if dm < 0:
            regrown_coords[:, 2] = -  regrown_coords[:, 2]

        regrown_coords = regrown_coords + np.array([a[0], a[1], a[2]])

        # if random.random() < .5:
        #     tmp = reverse(tmp)


        coordinates_list =  get_sequence_of_coords(self.length, i1, i2)
        cross_end  = coordinates_list[0]+ len(coordinates_list) -1 != coordinates_list[-1]

        if i1 < i2:
            if cross_end:
                # print('cross1')
                regrown_coords = regrown_coords.astype(float).T
            else:
                regrown_coords = regrown_coords.astype(float)[::-1].T
        else:
            if cross_end:
                # print('cross2')
                regrown_coords = regrown_coords.astype(float)[::-1].T
            else:
                regrown_coords = regrown_coords.astype(float).T



        self.coordinates[:, coordinates_list] = regrown_coords.copy()

        # print(i1, i2, self.length, coordinates_list, len(coordinates_list), regrown_coords.shape, a, b, n_beads)

        if  np.sum(np.abs(np.diff(self.coordinates.T, axis=0))) != len(np.diff(self.coordinates.T, axis=0)):
            print(np.sum(np.abs(np.diff(self.coordinates.T, axis=0))), len(np.diff(self.coordinates.T, axis=0)))
            print(repr(self.coordinates.T))
            print(repr(regrown_coords.T))
            sys.exit()

        self.output = self.coordinates.copy()
        # print(self.output)

        return self.output




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
        :type rotation:      Numpy array
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
        n_attempts = 0
        while  n_attempts<self.length:
            # selecting random bead
            n_attempts+=1

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
            if potentials:
                # selecting the axis with rotatios
                axis_rotation = random.choice(potentials)
                crank1 = random.choice(axis_rotation[0])
                rotation = random.choice(axis_rotation[1])  # index in the list of rotation matrices
                return crank, crank1, rotation

        return None, None, None



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
            if crank:
                self.coordinates = self.crankshaft_move(crank, crank1, rot[:, :, rotation])

        return self.coordinates


class Pin(Move):
    """
    class for  pin move
    """

    def __init__(self):
        super(Move, self).__init__()

    def __str__(self):
        return "pin"


    def make_move(self, n_steps):
        """
        runs pin move
        :return: (3,N) coords after performing pin move
        :rtype: numpy array,  int
        """

        # print('making pin move')
        pin, pin_plus, pin_minus = self.prepare_pin()

        # print("pin is", pin)
        if pin:
            # print('pin worked')
            self.coordinates = self.pin_move(pin, pin_plus, pin_minus)
        # else:
            # print('pin did not work')
        return self.coordinates



    def pin_move(self, pin, ping_plus, pin_minus):
        """
        performs pin move

        :return: (3,N) coords after performing pin move
        :rtype: numpy array,  int
        """

        new_position_delta = random.choice([[1,0,0], [-1,0,0], [0,1,0], [0,-1,0], [0,0,1], [0,0,-1]])

        self.coordinates[:, pin] = self.coordinates[:, pin_minus] + new_position_delta

        return  self.coordinates


    def prepare_pin(self):
        """
            prepares input data for the pin  move

        :return: 3-tuple as (preceding, actual, succeeding)
        :rtype: tuple of int
        """

        # the_same = False
        n_attempts = 0
        # pin_ = None

        while  (n_attempts < self.length):

            pin = np.random.randint(self.length)

            pin_minus = pin - 1
            pin_plus = pin + 1

            if pin == 0:
                pin_minus = self.length-1
            elif pin == self.length -1:
                pin_plus = 0

            if np.array_equal(self.coordinates[:, pin_minus], self.coordinates[:, pin_plus]):
                # the_same = True
                return pin,  pin_plus, pin_minus

            n_attempts +=1

        return   None, None, None




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
            if kink:
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

        n_attempts = 0

        while n_attempts < self.length:
            n_attempts += 1
            # print('make %i kinks' %n_attempts)
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

            if np.dot(v1, v2) == 0:
                return kink, kink_minus, kink_plus

        # print('%i kinks'% self.length)
        return None, None, None


