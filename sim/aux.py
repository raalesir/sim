"""
Auxiliary routines for  the modeler
-----------------------------------
"""
import  numpy as  np
import  random
import  math


try:
    from sim.consts import N, A, B, C
except ModuleNotFoundError:
    try:
        from consts import N, A,  B, C
    except:
        from .consts import N, A, B, C



def kink_move(coords, kink, kink_minus, kink_plus):
    """
            performs a kink move

    :param coords: 3D coordinates
    :type coords: numpy array
    :param kink: the bead subjected to kink move
    :type kink: int
    :param kink_minus: preceding bead
    :type kink_minus: int
    :param kink_plus: succeeding bead
    :type kink_plus: int
    :return: (3,N) coordinates after the kink move
    :rtype: Numpy array
    """


    if coords[0, kink_minus] == coords[0, kink] == coords[0, kink_plus]:
        coords[1, kink] = coords[1, kink_plus]
        coords[2, kink] = coords[2, kink_minus]
    # print('x')

    elif coords[1, kink_minus] == coords[1, kink] == coords[1, kink_plus]:
        coords[0, kink] = coords[0, kink_plus]
        coords[2, kink] = coords[2, kink_minus]
    # print('y')


    else:
        coords[0, kink] = coords[0, kink_plus]
        coords[1, kink] = coords[1, kink_minus]
    # print('z')

    return coords



def prepare_kink(coords_):
    """
        prepares input data for the kink  move

    :param coords_: 3D coordinates
    :type coords_: (3,N) Numpy array, int
    :return: 3-tuple as (preceding, actual, succeeding)
    :rtype: tuple of int
    """

    orthogonal = False

    while not orthogonal:

        kink = np.random.randint(N)

        kink_minus = kink - 1
        kink_plus = kink + 1

        if kink == 0:
            kink_minus = N - 1
        elif kink == N - 1:
            kink_plus = 0

            #     print(kink, length, coords[:, kink])

        v1 = coords_[:, kink] - coords_[:, kink_minus]
        v2 = coords_[:, kink_plus] - coords_[:, kink]

        orthogonal = np.dot(v1, v2) == 0

    return kink, kink_minus, kink_plus



def prepare_crank(coords_):
    """
       generates input data for crankshaft move --  crank points and the index in the rotation list

    :param coords_:    (3,N) coordinates
    :type coords_:  numpy array
    :return:    tuple of 3 integers, ``(crank, crank1, rotation)``, where ``crank`` is the starting and ``crank1`` the ending \
    indexes of the coordinates array for the crankshaft move: ``[crank+1:crank1-1]`` are subject to crankshaft; ``rotation`` -- \
    is the index of the list with rotation matrices
    :rtype:  tuple
    """

    potentials = [] # if this list is empty, it means that it the  chosen bead `crank` has to buddy
    # to form the rotation axis
    while not potentials:
        # selecting random bead
        crank = np.random.randint(N)
    #     print(crank)

        # generating potential positions for  the second bead
        potentials_x = np.where((coords_[1] == coords_[1, crank])  &  (coords_[2] == coords_[2, crank]))[0]
        potentials_y = np.where((coords_[0] == coords_[0, crank])  &  (coords_[2] == coords_[2, crank]))[0]
        potentials_z = np.where((coords_[0] == coords_[0, crank])  &  (coords_[1] == coords_[1, crank]))[0]


        potentials = [el for el in [(potentials_x, [0,1]),
                                    (potentials_y, [2,3]),
                                    (potentials_z, [4,5])] if len(el[0])>0
                     ]

        # selecting the axis with rotatios
        axis_rotation = random.choice(potentials)
        crank1 = random.choice(axis_rotation[0])
        rotation = random.choice(axis_rotation[1]) # index in the list of rotation matrices

    return crank, crank1, rotation





def crankshaft_move(coords_, A, B, rotation):
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

    tail = coords_[:, A+1:B]
    tail0 = coords_[:,A:A+1]
#     print(tail, tail0)
    tail = tail- tail0
#     print(tail)

#     print(rotation.shape)
    for i in range(tail.shape[1]):

        coords_[:, A+i+1] = np.matmul(rotation, tail[:,i]) + tail0[:, 0]

    return coords_



def make_circular_indexes(N):
    """

    :param N: number of bead
    :type N: int
    :return: tuple of 4 lists, each list keeps the indexes of beads belonging to a particular group
    :rtype: tuple
    """


    tmp = [i for i in range(N)]
    evens = tmp[::2]
    odds = tmp[1::2]
    tmp = len(evens)

    idx = int(np.ceil(N / 4))
    # print(idx)
    i0 = odds[-idx:]
    i1 = evens[:idx]
    i2 = evens[-tmp + idx:]
    i3 = odds[:tmp - idx]

    return i0, i1, i2, i3



def settle_init_point():
    """
    returns the initial point for the chain as the center of  the parallepiped specified by ``A`` and ``B``.


    :return: initial position of the polymer,  (1,3) numpy array
    :rtype: int
    """
    return np.array([A//2, B//2, C//2])



def make_circular_chain(N):
    """
        making circular chain of N beads

    :param N: number of beads
    :type N: int
    :return: *squashed* 3D-coordinates of the circular polymer
    :rtype: numpy array (3,N)
    """

    c = np.zeros((3, N))
    init_point = settle_init_point()

    i0, i1, i2, i3 = make_circular_indexes(N)
    print(i0, i1, i2, i3)
    c[:, i0] = (np.array([0, 1, 0]) + init_point).reshape(3, 1)
    c[:, i1] = (np.array([1, 1, 0]) + init_point).reshape(3, 1)
    c[:, i2] = (np.array([0, 0, 0]) + init_point).reshape(3, 1)
    c[:, i3] = (np.array([1, 0, 0]) + init_point).reshape(3, 1)

    return c



def generate_x_rotation_seq(N):
    """
        returns a sequence of multiples of 90 degrees  rotation around x-axis for a conformation \
        generated by ``make_circular_chain``. For example for ``N=12`` the output is ``[1, 1, 2, 1]``

    :param N: number  of beads
    :type N: int
    :return: ``[1,1,2,1,2,1,2,2,...]``
    :rtype: list
    """

    n_2 = 1
    n = 0
    seq = [[1]]

    while n < N / 2 - 2:
        t1 = [1] + n_2 * [2] + [1]
        t2 = t1[1:-1]
        n_2 += 1
        n = n + len(t1) + len(t2)
        seq.append(t1 + t2)

    return [item for sublist in seq for item in sublist][:int(N / 2 - 2)]


def make_rotation_matrices():
    """
        prepares rotation matrices

    :return: returns 6 matrices with rotations around X,Y,Z axis on plus/minus 90 degrees.
    :rtype: (3,3,6) numpy array
    """

    rot_z_1 = np.array([0, -1, 0, 1, 0, 0, 0, 0, 1]).reshape(3, 3, 1)
    rot_z_2 = np.array([0, 1, 0, -1, 0, 0, 0, 0, 1]).reshape(3, 3, 1)

    rot_x_1 = np.array([1, 0, 0, 0, 0, -1, 0, 1, 0]).reshape(3, 3, 1)
    rot_x_2 = np.array([1, 0, 0, 0, 0, 1, 0, -1, 0]).reshape(3, 3, 1)

    rot_y_1 = np.array([0, 0, 1, 0, 1, 0, -1, 0, 0]).reshape(3, 3, 1)
    rot_y_2 = np.array([0, 0, -1, 0, 1, 0, 1, 0, 0]).reshape(3, 3, 1)

    rot = np.concatenate([rot_x_1, rot_x_2, rot_y_1, rot_y_2, rot_z_1, rot_z_2], axis=2)

    return rot



def unroll_move(coords_, A, B, rotation):
    """
    unrolls the squashed conformation into a spiral one
    """

    # rotation = rot[:, :, 1]  # rot[:, :, np.random.randint(0, rot.shape[2])]

    tail = coords_[:, A + 1:B]
    tail0 = coords_[:, A:A + 1]
    tail = tail - tail0

    for i in range(tail.shape[1]):
        coords_[:, A + i + 1] = np.matmul(rotation, tail[:, i]) + tail0[:, 0]

    return coords_



def unroll_chain(coords, rotation_sequence, rotation):
    """
    unroll the  chain
    """

    for i in range(len(rotation_sequence)):
        for j in range(rotation_sequence[i]):
    #         print(i+1, N-i-2, j)
            coords = unroll_move(coords, i+1, N-i-2, rotation)
            # scatter.send_state('x'); scatter.send_state('y'); scatter.send_state('z')
            # lines.send_state('x');lines.send_state('y');lines.send_state('z')
    #         sleep(0.2)

    return coords



def run_kinks(coords, n_steps):
    """
    runs the ``n_steps`` random kinks
    :param coords:  (3,N) coords for ring polymer
    :type numpy array
    :param n_steps: number of kinks
    :type n_steps: int
    :return: (3,N) coords after performing ``n_steps`` random kinks
    :rtype:   numpy array, int
    """
    for i in range(n_steps):
        kink, kink_minus, kink_plus = prepare_kink(coords)
        coords = kink_move(coords, kink, kink_minus, kink_plus)

    return coords



def run_crankshafts(c, n_steps, rot):
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

    for i in range(n_steps):

        crank, crank1, rotation = prepare_crank(c)
        c = crankshaft_move(c, crank, crank1, rot[:, :, rotation])

    return c




def check_borders(c):
    """
        given coordinates ``c`` checks if any of the coordinates leaves the bound box

    :param c: 3D coordinates, (3,N)
    :type c: numpy array, int
    :return: True of False
    :rtype: bool
    """

    if (np.max(c[0,:]) < A) & (np.max(c[1,:]) < B) & (np.max(c[2]) < C) &\
        (np.min(c[0,:]) > 0) & (np.min(c[1,:]) >0) & (np.min(c[2]) >0):
        return True
    else:
        return False


def add(a,b):
    return a+b



def distance_noise_effects(d_, min_d):
    """
        learn how noise in the distance matrix affects the reconstructed distances

    :param d_: square matrix with distances squared
    :type d_: numpy 2D  array, float
    :param min_d: minimal distance used to filter the  matrix
    :type min_d: float
    :return:
    :rtype:
    """

    # difference with increasing the size of noise
    mu = 0.0
    delta = []
    ress = []
    N = d_.shape[1]
    #     sigmas = np.arange(0, 3.0, 0.05)
    sigmas = np.linspace(0, min_d / 2.0, 20)

    for sigma in sigmas:
        noise = np.random.normal(mu, sigma, size=(N, N))
        np.fill_diagonal(noise, 0.0)
        for i in range(N):
            for j in range(noise.shape[0]):
                noise[i, j] = noise[j, i]
                #         print(noise)
        tmp = np.sqrt(d_) + noise
        tmp[tmp < 0] = 0
        d_noise = np.multiply(tmp, tmp)
        m = Reconstruction.make_m(d_noise)
        u, s, vh = np.linalg.svd(m, full_matrices=False)
        res = np.matmul(u, np.sqrt(np.diag(s)))[:, :3]
        ress.append(res)
        # d_recovered = dist(res.T, cut=False)
        # d_recovered = d_recovered + d_recovered.T
        # delta.append((np.sum(np.sqrt(d_)) - np.sum(np.sqrt(d_recovered))) / np.sum(np.sqrt(d_)) * 100)

    return delta, sigmas, ress, noise, s



def loss_function(R, *args):

        tmp = np.array(R).reshape((4, 3))
        rot = tmp[:3,:]
        v = tmp[-1, :]

        r0 = args[0]
        r1 = args[1] + v
        r1 = v + np.matmul(r1-v, rot)
        configuration_distance = sum([
                                         np.dot(tmp, tmp) for tmp in r0 - r1
                                         ])
        return np.sqrt(configuration_distance / r0.shape[0])



def make_m(d):
    """
    making Gram matrix from distance
    """
    N = d.shape[1]
    m = np.zeros((N, N))

    for i in range(N):
        for j in range(N):
            m[i, j] = (d[0, i] + d[j, 0] - d[i, j]) / 2.
    return m


def n_conf(N, dx, dy, dz):
    """
    calculates the number of conformations  of ideal grid polymer given
    number of bonds and displacements along the grid.
    """
    dx = abs(dx); dy = abs(dy); dz = abs(dz)


    if ((N - dx - dy + dz) % 2 != 0) | ((N - dx - dy - dz) % 2 != 0):
        return 0
    else:

        n_plus = int((N - dx - dy + dz) / 2)
        n_minus = int((N - dx - dy - dz) / 2)

        numerator = math.factorial(N)
        res = 0.0
        for x in range(n_minus + 1):
            for y in range(n_minus - x + 1):
                res += numerator / math.factorial(x) / math.factorial(x + dx) / math.factorial(y) / math.factorial(
                    y + dy) / \
                       math.factorial(n_plus - x - y) / math.factorial(n_minus - x - y)

        return res



def cache_n_conf(N, dx, dy, dz):
    """
    caches the n_conf for each point on the grid given by dz, dy, dz
    """
    res = []
    for n in range(N):
        for i in range(dx):
            for j in range(dy):
                for k in range(dz):
                    res.append(n_conf(n + 1, i, j, k))

    return np.array(res).reshape(N, dx, dy, dz)  # .astype(int)



def get_n_beads(n, i1 ,i2):
    """
    given length of the ring polymer (n) and the starting-ending positions (i1, i2)
    returns the length of the segment of the smaller length

    :param n: number of elements
    :type n: int
    :param i1: starting position
    :type i1: int
    :param i2: ending position
    :type i2: int
    :return: length of the smaller segment
    :rtype: int
    """
    if (n < max(i1,i2)) | (i1<0) | (i2<0):
        return  0
    else:
        return min(abs(i2-i1)+1, n - abs(i2-i1)+1)



def get_sequence_of_coords1(n, i1, i2):
    """
    given  number  of elements in the ring (n), starting and ending indexes (i1, i2)
    returns a list of element numbers subjected to move

    :param n: number of elements
    :type n: int
    :param i1: starting position
    :type i1: int
    :param i2: ending position
    :type i2: int

    :return: list of elements of ring chain subjected to move
    :rtype: list
    """

    if (n < max(i1, i2)) | (i1 < 0) | (i2 < 0):
        return None
    else:
        if abs(i2-i1) < n - abs(i2-i1):
            if i1<i2: return  list(range(i1, i2+1))
            else: return  list(range(i2, i1+1))
        else:
            if i1<i2: return list(range(i2, n)) + list(range(i1+1))
            else: return list(range(i1, n)) + list(range(i2+1))



def get_sequence_of_coords(n, i1, i2, ori_ark = True):
    """
    given  number  of elements in the ring (n), starting and ending indexes (i1, i2)
    returns a list of element numbers subjected to move

    :param n: number of elements
    :type n: int
    :param i1: starting position
    :type i1: int
    :param i2: ending position
    :type i2: int
    :param ori_ark: return ark passing ORI
    :type ori_ark: bool
    :return: list of elements of ring chain subjected to move
    :rtype: list
    """

    if (n < max(i1, i2)) | (i1 < 0) | (i2 < 0):
        return None
    else:

        if i1<i2:
            r1 = range(i1,i2+1)
            r2 = range(i2, n)
            r3 = range(i1+1)
        else:
            r1 = range(i2, i1 + 1)
            r2 = range(i1, n)
            r3 = range(i2 + 1)


        # if len(r1) <= len(r2) +  len(r3):
        #         small = list(r1); large = list(r2) + list(r3)
        # else:
        #         large = list(r1)
        #         small = list(r2) + list(r3)


        if ori_ark:
            return  list(r2) + list(r3)
        else:
            return  list(r1)
