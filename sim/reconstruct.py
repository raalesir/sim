"""
    Reconstruction of  Coarse-Grained structure  of the polymer
    ==============================================================


"""
import  numpy as  np
import  random
from scipy.optimize import minimize
import  matplotlib.pyplot as plt
import  sys

from sim.sim.aux import loss_function, make_m


class Reconstruction:
    """
    Takes  the  3D coordinates of ring grid polymer and performs  analysis of the conformation
    """

    def __init__(self, coords):
        """

        :param coords: 3D coords on cubic grid
        :type coords: Numpy array, int
        """
        # np.set_printoptions(threshold=sys.maxsize)

        self.coords = coords
        self.distance_matrix = self.get_dist() # distance matrix


    def  __str__(self):
        r  = 40*"X" + """
         coordinates shape: %s
"""%str(self.coords.shape) + 40*"X"

        return  r


    def get_dist(self):
        """
        squared distances from coordinates

        """
        n = self.coords.shape[1]
        d = np.zeros((n, n))
        for i in range(n):
            for j in range(i):
                tmp = self.coords[:, i] - self.coords[:, j]
                d[i, j] = np.dot(tmp, tmp)
        # if cut:
        #     # print('cutoff %i' % cutoff)
        #     d[d < cutoff] = 0

        return  d


    def filter_distance_matrix(self, cutoff=100):
        """
        filters the distance matrix if needed

        :return: self.distance_matrix_filtered
        :rtype: Numpy array, float
        """
        self.distance_matrix[self.distance_matrix <  cutoff] = 0



    def create_coarse_coords(self, cutoff_list):
        """
        for each of threshould in the `cutoff_list`  finds all  possible sets of `self.coords` with max
        number  of beads

        :param cutoff_list: set of cutoff for the distance matrix
        :type cutoff_list: list of floats
        :return: dict with keys as elements from `cutoff_list` and the values as list of Numpy arrays with
        coordinates  satisfying the cutoff.
        :rtype: dict
        """

        coarse_coords = {}
        print('cutoff  | number of coords in the set | # of points in coords')
        for cutoff_ in cutoff_list:

            d = np.where(self.distance_matrix >= cutoff_, self.distance_matrix, 0)
            top, bottom = np.nonzero(d)

            el_top = sorted(np.unique(top))[::-1]

            coarse_indexes = self.get_max_coarse_grained_beads(el_top, top, bottom)
            print(cutoff_, len(coarse_indexes), len(coarse_indexes[0]))

            tmp = []
            for subitem in coarse_indexes:
                tmp.append(self.coords[:, subitem])
            coarse_coords[cutoff_] = tmp

        return coarse_coords


    @staticmethod
    def dist(c):
        """
        squared distances from coordinates

        """
        n = c.shape[1]
        d = np.zeros((n, n))
        for i in range(n):
            for j in range(i):
                tmp = c[:, i] - c[:, j]
                d[i, j] = np.dot(tmp, tmp)

        return d


    def create_rsme(self,  coarse_coords):
        rmse = {}
        for cutoff_ in coarse_coords:
            print('the cutoff is %f; M=%i' % (cutoff_, coarse_coords[cutoff_][0].shape[1]))
            # print(coarse_coords[cutoff_][0].shape)
            d = Reconstruction.dist(coarse_coords[cutoff_][0])
            # print(d)
            # d = np.where(self.distance_matrix >= cutoff_, self.distance_matrix, 0)
            d = d + d.T

            tmp_precision = []
            tmp_processed_coords =  []
            c_off = np.sqrt(cutoff_)
            for series in range(10):
                delta, sigmas, ress, noise, ss = Reconstruction.distance_noise_effects(d, c_off)
                # print('ress', ress)
                processed_coords, precision = Reconstruction.align_results(ress)
                tmp_precision.append(precision)
                tmp_processed_coords.append(processed_coords)
            rmse['d_min= ' + "{:.2f}".format(c_off) + ' ;__M=' + str(coarse_coords[cutoff_][0].shape[1])] = \
                np.array(tmp_precision).mean(axis=0), np.array(tmp_processed_coords).mean(axis=0)


        return rmse



    @staticmethod
    def align_results(raw_results):
        """

        :param raw_results: list  of restored coordinates of the same config as a function of `sigma`
        :type raw_results: list
        :return: tuple of lists with aligned coords and corresponding precision
        :rtype: tuple
        """
        processed = []
        precision = [0]
        processed.append(raw_results[0])
        for i in range(1, len(raw_results)):
            result = minimize(loss_function, (0, 0, 0, 0, 0, 0, 0, 0, 0,0,0,0), args=(raw_results[0], raw_results[i]),
                              method='CG')
            #         print(result.fun, loss_function(np.diag((1,1,1)), raw_results[0], raw_results[i]))
            precision.append(result.fun)
            res = result.x.reshape((4, 3))
            rot  = res[:3,:] # rotation matrix
            v = res[-1, :] # parallel transfer
            tmp = np.matmul(raw_results[i], rot) + v

            processed.append(tmp)
        return processed, precision


    @staticmethod
    def get_sigmas(x, y, A, M,  n_T=1):
        """
        returns periodic noise
        with given amplitude
        """

        xx, yy = np.meshgrid(x, y)
        res = np.abs(
            # np.random.normal(0, A/10 , size=(M, M)) + \
              A * np.sin(n_T * 2 * np.pi * xx / M) * np.sin(n_T * 2 * np.pi * yy / M)
        )#+ A
        res[np.abs(res) < 0.5*A] = 0

        return res


    @staticmethod
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
        sigmas = np.linspace(0, min_d / 2.0, 30)

        for sigma in sigmas:
            # print('sigmma', sigma)
            # noise = np.random.normal(mu, sigma, size=(N, N))
            noise = np.random.normal(mu, Reconstruction.get_sigmas(
                np.arange(1,N+1), np.arange(1,N+1), sigma, N
            ), size=(N, N))

            # noise = Reconstruction.get_noise(np.arange(1,N+1), np.arange(1,N+1), sigma*2, N)
            np.fill_diagonal(noise, 0.0)
            for i in range(N):
                for j in range(noise.shape[0]):
                    noise[i, j] = noise[j, i]
                    #         print(noise)
            tmp = np.sqrt(d_) + noise
            tmp[tmp < 0] = 0
            d_noise = np.multiply(tmp, tmp)
            m = make_m(d_noise)
            u, s, vh = np.linalg.svd(m, full_matrices=False)
            res = np.matmul(u, np.sqrt(np.diag(s)))[:, :3] # dim  = (M,3)  instead of (3,M)
            ress.append(res)
            # d_recovered = dist(res.T, cut=False)
            # d_recovered = d_recovered + d_recovered.T
            # delta.append((np.sum(np.sqrt(d_)) - np.sum(np.sqrt(d_recovered))) / np.sum(np.sqrt(d_)) * 100)

        return delta, sigmas, ress, noise, s


    @staticmethod
    def unroll_indexes(idx):
        """
        unrolls  [192, [172, [94, [70]]]]  to [192, 172, 94, 70] """

        go = True
        unrolled_indexes = []
        while go:
            unrolled_indexes.append(idx[0])
            if len(idx) > 1:
                idx = idx[1]
            else:
                go = False
        return unrolled_indexes


    @staticmethod
    def get_max_coarse_grained_beads(el_top, top, bottom):
        """
        returns the first  set (if many) of coordinate  indexes having the minimal mutual distance
        """


        def find_beads2(current_node, search_space, found):
            """
            finds  the child in the `search_space` for a `current_node` with the most number of  children.
            Returns the first one if there  are several children with equal number of descedents
            """
            # https://stackoverflow.com/questions/39604394/counting-paths-in-a-recursive-call
            #     print("%s --> %s; %s"%(max_, search_space))

            if len(search_space) < 1:
                found[current_node] = 0
                return [current_node]

            elif current_node in found:
                return found[current_node]

            else:
                node_most_children = 0
                most_num_children = []  # 0
                tmp = 0
                for child in search_space[::-1]:
                    search_space_child = bottom[np.where(top == child)]
                    search_space = np.intersect1d(search_space, search_space_child)

                    tmp = [current_node] + [find_beads2(child, search_space, found)]

                    if len(tmp) > len(most_num_children):
                        most_num_children = tmp
                        node_most_children = child
                        found[current_node] = tmp

                return most_num_children

        collect_deepest = []
        max_depth = []  # 0
        for max_ in el_top[:]:
            node_most_child = []
            depth = find_beads2(max_, bottom[np.where(top == max_)], {})
            #         print(len(node_most_child), depth)
            #         print(node_most_child)

            ids = Reconstruction.unroll_indexes(depth)
            if len(ids) >= len(max_depth):
                max_depth = ids
                collect_deepest.append(ids)

        deepest = [el for el in collect_deepest if len(el) >= len(max_depth)]

        return deepest #random.choice(deepest)


    @staticmethod
    def plot_rmse(rmse):

        plt.figure(figsize=(16, 12))

        for key in rmse:
            min_dist = key.split(' ')[1]
            #     sigmas_ = [el/float(min_dist) for el in sigmas]
            t = [el / float(min_dist) for el in rmse[key][0]]
            sigmas = np.linspace(0, 0.5, 30)

            plt.plot(sigmas, t, '.-',linewidth=5, label=key)

        plt.title("Relative deviation of disturbed from reconstructed  structures", fontsize=20)
        plt.xlabel(r"$\sigma/d_{min}$", fontsize=20)
        plt.ylabel("$RMSE/d_{min}$ from the original structure", fontsize=20)
        plt.grid()
        # plt.xlim(0,0.4)
        # plt.ylim(0, 1)
        # plt.xscale('log')

        plt.legend()
        plt.savefig('loss1.png')
        plt.show()



if __name__ ==  "__main__":

    c = np.array([[0,0,1,1], [1,0,0,1], [0,0,0,0]])
    c  = np.array([[16., 17., 17., 16., 15., 15., 16., 17., 18., 18., 18., 19., 19.,
        19., 20., 20., 20., 20., 20., 20., 19., 19., 19., 19., 20., 20.,
        19., 19., 19., 19., 18., 18., 18., 18., 18., 18., 17., 17., 17.,
        16., 16., 15., 15., 15., 14., 14., 13., 13., 13., 12., 11., 10.,
        10., 10., 10., 10., 10., 10., 10., 11., 11., 11., 12., 12., 13.,
        14., 14., 15., 16., 16., 17., 17., 18., 18., 18., 18., 18., 17.,
        16., 16., 17., 17., 17., 17., 17., 17., 17., 17., 16., 16., 15.,
        14., 14., 15., 15., 15., 14., 14., 14., 14., 14., 13., 13., 13.,
        13., 12., 12., 12., 13., 14., 14., 15., 15., 16., 16., 16., 16.,
        15., 15., 15., 15., 15., 15., 14., 14., 14., 13., 12., 12., 13.,
        13., 14., 15., 15., 15., 15., 14., 14., 15., 15., 15., 15., 15.,
        14., 14., 14., 15., 15., 16., 16., 17., 17., 17., 17., 18., 18.,
        17., 17., 17., 16., 16., 16., 15., 15., 16., 16., 16., 15., 14.,
        13., 13., 14., 15., 15., 15., 15., 15., 15., 15., 15., 14., 13.,
        13., 13., 12., 12., 12., 12., 12., 12., 12., 13., 13., 13., 14.,
        14., 14., 14., 15., 16.],
       [16., 16., 15., 15., 15., 15., 15., 15., 15., 14., 14., 14., 15.,
        15., 15., 15., 16., 17., 17., 18., 18., 18., 18., 19., 19., 20.,
        20., 20., 19., 19., 19., 18., 18., 19., 20., 20., 20., 20., 20.,
        20., 20., 20., 21., 21., 21., 22., 22., 21., 21., 21., 21., 21.,
        21., 20., 20., 19., 18., 18., 19., 19., 19., 19., 19., 19., 19.,
        19., 19., 19., 19., 19., 19., 18., 18., 17., 16., 15., 15., 15.,
        15., 15., 15., 14., 13., 12., 12., 11., 10., 10., 10., 10., 10.,
        10.,  9.,  9.,  9.,  9.,  9.,  9.,  9., 10., 10., 10., 10.,  9.,
         8.,  8.,  8.,  8.,  8.,  8.,  9.,  9.,  8.,  8.,  9., 10., 11.,
        11., 12., 13., 13., 12., 12., 12., 13., 13., 13., 13., 12., 12.,
        12., 12., 12., 13., 14., 14., 14., 15., 15., 16., 16., 17., 17.,
        17., 17., 17., 17., 17., 17., 17., 17., 18., 18., 18., 18., 18.,
        18., 19., 20., 20., 19., 19., 19., 19., 19., 20., 20., 20., 20.,
        20., 21., 21., 21., 21., 21., 20., 20., 19., 19., 19., 19., 19.,
        18., 17., 17., 16., 16., 16., 16., 15., 14., 14., 15., 16., 16.,
        16., 16., 16., 16., 16.],
       [15., 15., 15., 15., 15., 14., 14., 14., 14., 14., 15., 15., 15.,
        16., 16., 17., 17., 17., 16., 16., 16., 17., 18., 18., 18., 18.,
        18., 19., 19., 20., 20., 20., 19., 19., 19., 20., 20., 19., 18.,
        18., 17., 17., 17., 16., 16., 16., 16., 16., 17., 17., 17., 17.,
        16., 16., 17., 17., 17., 18., 18., 18., 19., 20., 20., 21., 21.,
        21., 20., 20., 20., 21., 21., 21., 21., 21., 21., 21., 22., 22.,
        22., 21., 21., 21., 21., 21., 20., 20., 20., 21., 21., 20., 20.,
        20., 20., 20., 19., 18., 18., 17., 16., 16., 15., 15., 14., 14.,
        14., 14., 13., 12., 12., 12., 12., 12., 12., 12., 12., 12., 12.,
        12., 12., 12., 11., 11., 10., 10., 10.,  9.,  9.,  9.,  9.,  9.,
         8.,  8.,  8.,  8.,  8.,  7.,  7.,  7.,  7.,  7.,  6.,  6.,  5.,
         5.,  4.,  3.,  3.,  4.,  4.,  3.,  3.,  3.,  4.,  5.,  5.,  6.,
         6.,  6.,  6.,  6.,  6.,  5.,  5.,  4.,  4.,  4.,  3.,  3.,  3.,
         3.,  3.,  3.,  3.,  4.,  5.,  5.,  6.,  6.,  7.,  8.,  8.,  8.,
         8.,  8.,  8.,  8.,  9., 10., 11., 11., 11., 11., 11., 11., 11.,
        12., 13., 14., 14., 14.]])


    r = Reconstruction(c)
    print(r)
    print(r.distance_matrix)
    # r.filter_distance_matrix(cutoff=100)
    # print(r.distance_matrix)
    coarse_coordinates = r.create_coarse_coords(list(range(1, 100, 20))
                                                )
    rmse = r.create_rsme(coarse_coordinates)
    r.plot_rmse(rmse)



