"""
    Reconstruction of  Coarse-Grained structure  of the polymer
    ==============================================================


"""
import  numpy as  np
import  random


class Reconstruction:
    """
    Takes  the  3D coordinates of ring grid polymer and performs  analysis of the conformation
    """

    def __init__(self, coords):
        """

        :param coords: 3D coords on cubic grid
        :type coords: Numpy array, int
        """

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
                tmp = c[:, i] - c[:, j]
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
        for cutoff_ in cutoff_list:
            top, bottom = np.nonzero(self.distance_matrix)
            el_top = sorted(np.unique(top))[::-1]

            coarse_indexes = self.get_max_coarse_grained_beads(el_top, top, bottom)
            print(len(coarse_indexes))
            coarse_coords[cutoff_] = self.coords[:, coarse_indexes]

        return coarse_coords



    @staticmethod
    def get_max_coarse_grained_beads(el_top, top, bottom):
        """
        returns the first  set (if many) of coordinate  indexes having the minimal mutual distance
        """

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

            ids = unroll_indexes(depth)
            if len(ids) >= len(max_depth):
                max_depth = ids
                collect_deepest.append(ids)

        deepest = [el for el in collect_deepest if len(el) >= len(max_depth)]

        return random.choice(deepest)





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
    r.filter_distance_matrix(cutoff=100)
    print(r.distance_matrix)
    print(r.create_coarse_coords(range(1, 50, 4)))



