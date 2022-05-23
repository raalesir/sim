"""
    calculates number of configurations for a ring grid phantom polymer and
    overlap distribution for the chain.
"""

import math
from math import comb
from collections import Counter
import json
import os
import matplotlib.pyplot as plt


class Overlap:
    """
    calculating number of overlaps  for a ring grid polymer of given N.
    N should be  even to make a ring polymer.
    """

    def __init__(self, N):
        """
        N -- number of monomers
        """
        self.n = N
        self.n_conformations = self.n_conform()
        self.overlaps_hist = None
        self.indexes = self.calculate_steps()
        self.dict = self.make_steps()
        self.encoded_conformations = None

    #         self.keep_result = []




    def __str__(self):
        return "overlaps for %i beads has %i conformations" % (self.n, self.n_conformations)

    def n_conform(self):
        """
        """
        r = 0
        for i in range(self.n // 2 + 1):
            for j in range(self.n // 2 - i + 1):
                r = r + math.factorial(self.n) / math.factorial(i)**2 / math.factorial(j)**2 / math.factorial(
                    self.n // 2 - i - j)**2
        return r


    def fun(self, d, res):
        """
        returns a letter representation for a number combination.
        For example for  n=4 it  is (2,0,0) -->> ['i+i+i-i-', 'i+i-i+i-', 'i+i-i-i+', '-+-+', '-++-'  ,'--++']

        :param d: numeric representation of the closed trajectory e.g. (2,0,1)
        :type d: dict
        :param res: string representation for all possible configurations  for the given numeric representation, see above in the function doc
        :type res: str
        :return: string representation for all possible configurations
        :rtype:
        """
        if sum(d.values()) == 0: # if (0,0,0) that is we are at the initial point
            Overlap.keep_result(res)
            return
        else:
            for k in [item for item in d.keys() if d[item] > 0]:
                r = res
                r += k

                tmp = d.copy()
                tmp[k] -= 1

                self.fun(tmp, r)


    def calculate_steps(self):
        """
        given number of monomers, n, produce the indexes (i,j,k)
        as the number of steps to make in positive and negative direction
        """
        res = []
        for i in range(self.n // 2 + 1):
            for j in range(self.n // 2 - i + 1):
                res.append((i, j, self.n // 2 - i - j))
        return res


    def make_step(self, tup):
        """
        encodes  single index
        :return:
        :rtype:
        """
        res = []
        d = {}
        d['i+'] = tup[0]
        d['i-'] = tup[0]
        d['j+'] = tup[1]
        d['j-'] = tup[1]
        d['k+'] = tup[2]
        d['k-'] = tup[2]

        res.append(d)

        return res


    def make_steps(self):
        """
        encodes indexes in dict
        """
        res = []
        for tup in self.indexes:
            d = {}
            d['i+'] = tup[0]
            d['i-'] = tup[0]
            d['j+'] = tup[1]
            d['j-'] = tup[1]
            d['k+'] = tup[2]
            d['k-'] = tup[2]

            res.append(d)
        return res

    #     @static
    def keep_result(data):
        """
        """
        Overlap.keep_result.all.append(data)


    def calculate_all_conformations(self):
        Overlap.keep_result.all = []

        for entry in self.dict:
            # generating trajectory representation  for each combination
            self.fun(entry, '')
            # print(Overlap.keep_result.all)


    def encode_single_conformation(self, conformation):
        """

        :param conformation:
        :type conformation:
        :return:
        :rtype:
        """

        conf_encoded = []
        start = [0, 0, 0]
        for symbol in [conformation[i:i + 2] for i in range(0, len(conformation), 2)]:
            if symbol == 'k+':
                start[2] += 1
            elif symbol == 'k-':
                start[2] -= 1
            elif symbol == 'i+':
                start[0] += 1
            elif symbol == 'i-':
                start[0] -= 1
            elif symbol == 'j+':
                start[1] += 1
            elif symbol == 'j-':
                start[1] -= 1
            conf_encoded.append(tuple(start.copy()))

        return conf_encoded


    def encode_to_coords(self):
        """
        encodes string conformation representation to the numeric representation
        e.g. 'i+i+i-i-' --> [ (1,0,0), (2,0,0), (1,0,0), (0,0,0) ]
        """
        res = []
        for conformation in Overlap.keep_result.all:
            conf_encoded = []
            start = [0, 0, 0]
            for symbol in [conformation[i:i + 2] for i in range(0, len(conformation), 2)]:
                if symbol == 'k+':
                    start[2] += 1
                elif symbol == 'k-':
                    start[2] -= 1
                elif symbol == 'i+':
                    start[0] += 1
                elif symbol == 'i-':
                    start[0] -= 1
                elif symbol == 'j+':
                    start[1] += 1
                elif symbol == 'j-':
                    start[1] -= 1
                conf_encoded.append(tuple(start.copy()))

            res.append(conf_encoded)
        self.encoded_conformations = res




    def get_overlaps(self):
        """
        """
        overlaps = []
        for conf in self.encoded_conformations:
            #  conf is like [ (1,0,0), (2,0,0), (1,0,0), (0,0,0) ]
            overlaps.append(sum([comb(lst, 2) for lst in Counter(conf).values()]))

        counts = Counter(overlaps)

        self.overlaps_hist = dict(counts)


    def get_overlaps_histogram(self):

        fname = "counts_%i.json" % (self.n)
        if not os.path.isfile(fname):

            self.calculate_all_conformations() # all confs are encoded as a list of strings like 'i+i+i-i-'
            self.encode_to_coords()
            self.get_overlaps()

        else:
            dct = open(fname, 'r').read()
            dct = json.loads(dct)
            self.overlaps_hist = dict(zip([int(el) for el in dct.keys()], dct.values()))

        return self.overlaps_hist


    def save_overlaps_histogram(self):
        fname = "counts_%i.json" % (self.n)
        if not os.path.isfile(fname):
            json.dump(self.overlaps_hist, open(fname, 'w'))


    def plot_overlaps_histogram(self):

        self.overlaps_hist = dict(zip(self.overlaps_hist.keys(), [v/sum(self.overlaps_hist.values()) for v in self.overlaps_hist.values()]))
        plt.bar(self.overlaps_hist.keys(), self.overlaps_hist.values())
        plt.yscale('log')
        plt.xlabel('number of overlaps')
        plt.ylabel('number of conformations')



if __name__=="__main__":

    overlaps = Overlap(8)
    print(overlaps)

    print(overlaps.get_overlaps_histogram())

    overlaps.save_overlaps_histogram()

    overlaps.plot_overlaps_histogram()
