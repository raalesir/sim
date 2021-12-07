"""
    Polymer class
    =============

"""

import  numpy as np


try:
    from sim.cell import CubicCell
    from sim.moves import  Kink, CrankShaft, Pin
    from sim.consts import  rot
    from sim.aux import get_sequence_of_coords
except ModuleNotFoundError:
    try:
        from cell import CubicCell
        from moves import Kink, CrankShaft, Pin
        from consts import rot
        from aux import get_sequence_of_coords
    except:
        from .cell import CubicCell
        from .moves import Kink, CrankShaft, Pin
        from .consts import rot
        from  .aux import get_sequence_of_coords
# except:
#     from  sim.cell import  CubicCell



class Polymer:
    """
    represents a ring polymer on a cubic grid
    """

    def __init__(self, n,  cell_, *args, **kwargs):
        """
        polymer init
        """
        self.n = n
        self.cell = cell_
        self.md_size = int(self.n // 4)

        self.coords_tmp = None
        self.coords = self.set_init_coords()

        for move in kwargs.pop('moves'):
            name = 'move_' + move.__str__()
            setattr(self, name, move)

        for md in kwargs.pop('mds'):
            if md == 'ori':
                self.ori_md = self.get_ori_md()

        # self.kink.coordinates = self.coords
        # self.kink.length = self.n

        self.distance_matrix = None



    def __str__(self):
        msg = 40*"X" + """\nring polymer with %i beads\nconfined into  cell with the dimensions (%i,%i,%i)\nmoves registered are: %s\n"""\
               %(self.n, self.cell.A, self.cell.B, self.cell.C, [move  for  move in dir(self) if 'move_' in move]) +\
              40*"X"

        return msg



    def get_ori_md(self):
        """
        returns list of indexes for ORI Macro Domain based on the chain length

        :return: [A, B] first and last bead of the ORI macrodomain
        :rtype: list
        """

        md_size = self.md_size # approximate size
        i1 = self.n - md_size//2
        i2 = md_size//2
        coordinates_list = get_sequence_of_coords(self.n, i1, i2, ori_ark=True)

        return coordinates_list


    def get_cm_md(self, tmp=False):
        """
        returns center of mass for the Macro Domain

        :return:
        :rtype:
        """
        if tmp:
            md = self.coords_tmp[:, self.ori_md]
        else:
            md = self.coords[:, self.ori_md]

        return np.mean(md, axis=1)



    def check_borders(self):
        """
            given coordinates ``c`` checks if any of the coordinates leaves the bound box

        :return: True of False
        :rtype: bool
        """

        if (np.max(self.coords_tmp[0, :]) < self.cell.A) &\
                (np.max(self.coords_tmp[1, :]) < self.cell.B) &\
                (np.max(self.coords_tmp[2]) < self.cell.C) & \
                (np.min(self.coords_tmp[0, :]) >= 0) & (np.min(self.coords_tmp[1, :]) >= 0) &\
                (np.min(self.coords_tmp[2]) >= 0):
            return True
        else:
            return False


    def check_overlap(self):
        """
        checking chain for overlaps

        :return: True/False
        :rtype: bool
        """

        return self.coords_tmp.shape == np.unique(self.coords_tmp, axis=1).shape


    def calculate_distances(self):
        """
        calculating distance matrix

        :return: distance matrix `NxN`
        :rtype: numpy array
        """

        size_ = self.coords.shape[1]
        dst = np.zeros((size_, size_))
        for i in range(size_):
            for j in range(0, i):
                tmp = self.coords[:, i] - self.coords[:, j]
                dst[i, j] = np.dot(tmp, tmp)

        self.distance_matrix = dst

    # def save_distan

    def _make_circular_indexes(self):
        """

            :param N: number of bead
            :type N: int
            :return: tuple of 4 lists, each list keeps the indexes of beads belonging to a particular group
            :rtype: tuple
            """

        tmp = [i for i in range(self.n)]
        evens = tmp[::2]
        odds = tmp[1::2]
        tmp = len(evens)

        idx = int(np.ceil(self.n / 4))
        # print(idx)
        i0 = odds[-idx:]
        i1 = evens[:idx]
        i2 = evens[-tmp + idx:]
        i3 = odds[:tmp - idx]

        return i0, i1, i2, i3


    def set_init_coords(self):
        """
            making circular chain of N beads

       :return: *squashed* 3D-coordinates of the circular polymer
       :rtype: numpy array (3,N)
        """

        c = np.zeros((3, self.n))
        init_point = self.cell.get_center()

        i0, i1, i2, i3 = self._make_circular_indexes()

        # print(i0, i1, i2, i3)
        c[:, i0] = (np.array([0, 1, 0]) + init_point).reshape(3, 1)
        c[:, i1] = (np.array([1, 1, 0]) + init_point).reshape(3, 1)
        c[:, i2] = (np.array([0, 0, 0]) + init_point).reshape(3, 1)
        c[:, i3] = (np.array([1, 0, 0]) + init_point).reshape(3, 1)

        return c



    def get_coords(self):
        return self.coords


    def generate_x_rotation_seq(self):
        """
            returns a sequence of multiples of 90 degrees  rotation around x-axis for a conformation \
            generated by ``make_circular_chain``. For example for ``N=12`` the output is ``[1, 1, 2, 1]``


        :return: ``[1,1,2,1,2,1,2,2,...]``
        :rtype: list
        """

        n_2 = 1
        n = 0
        seq = [[1]]

        while n < self.n / 2 - 2:
            t1 = [1] + n_2 * [2] + [1]
            t2 = t1[1:-1]
            n_2 += 1
            n = n + len(t1) + len(t2)
            seq.append(t1 + t2)

        return [item for sublist in seq for item in sublist][:int(self.n / 2 - 2)]



    def unroll_move(self, A, B, rotation):
        """
        unrolls the squashed conformation into a spiral one
        """

        # rotation = rot[:, :, 1]  # rot[:, :, np.random.randint(0, rot.shape[2])]

        tail = self.coords[:, A + 1:B]
        tail0 = self.coords[:, A:A + 1]
        tail = tail - tail0

        for i in range(tail.shape[1]):
            self.coords[:, A + i + 1] = np.matmul(rotation, tail[:, i]) + tail0[:, 0]

        return self.coords


    def unroll_chain(self):
        """
        unroll the  chain
        """

        rotation_sequence = self.generate_x_rotation_seq()
        rotation = rot[:, :, 1] #self._make_rotation_matrices()[:, :, 1]
        # print("rotation_sequence", rotation_sequence)

        for i in range(len(rotation_sequence)):
            for j in range(rotation_sequence[i]):
                self.coords = self.unroll_move(i + 1, self.n - i - 2, rotation)

        return self.coords



if __name__ == "__main__":


    cell = CubicCell(6,10,6)
    print(cell)
    print(cell.get_center())

    kink = Kink()
    print(kink)
    crankshaft = CrankShaft()
    print(crankshaft)

    pin = Pin()
    print(pin)


    polymer = Polymer(20, cell, kink, crankshaft, pin)
    print(polymer)
    print('squashed coordinates are:')
    print(polymer.get_coords())

    print('unrolled coordinates are:')
    print(polymer.unroll_chain())
    # print('polymer coords\n', polymer.coords)

    if hasattr(polymer, 'move_kink'):
        print('kinking...')
        kink.coordinates = polymer.coords#.copy()
        polymer.coords =  polymer.move_kink.getOutput(100)
        print(polymer.coords)

    if hasattr(polymer, 'move_crankshaft'):
        print('crankshafting...')
        crankshaft.coordinates = polymer.coords
        polymer.coords = polymer.move_crankshaft.getOutput(100)
        print(polymer.coords)






