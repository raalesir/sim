"""
    Class for a  force-field.
    The field  acts  on the chain causing it to move towards the "center of gravity",
    minimizing its  potential energy.

"""
import  numpy as np


class ForceField:

    def __init__(self, on=True, amplitude=1, linear=False):
        """

        """
        self.on = on # turn on/off the potential
        self.linear = linear

        self.origin = self.get_origin()
        self.f_r = self.get_value()
        self.amplitude = amplitude


    def __str__(self):
        if self.linear:
            return 'linear'
        else:
            return 'Euclidean'


    def get_value(self):
        """
        potential dependency on  distance, electrostatic-like

        :return: function for  potential
        :rtype: function
        """
        def f(x):
            return  -self.amplitude/(x+0.01)

        return   f


    def get_origin(self):
        if self.linear:
            return 0
        else:
            return  0,0,0



    def get_distance(self, a, b=None):
        """
        returns  distance

        :return:
        :rtype:
        """

        if self.linear:
            return   self.get_lin_distance(a, b)

        else:
            return  self.get_eu_distance(a, b)



    def get_eu_distance(self, a, b=None):
        """
        calculates Euclidean squared  distance between points

        :param a: point A= (k1,l1,m1)
        :type a: tuple
        :param b: point B= (k2,l2,m2)
        :type b: tuple
        :return: distance
        :rtype: float
        """
        if b:
            tmp = [a[0] - b[0], a[1]-b[1], a[2]-b[2]]
        else:
            o = self.origin
            tmp = [a[0] - o[0], a[1]-o[1], a[2]-o[2]]

        return  np.dot(tmp,tmp)



    def get_lin_distance(self, a, b = None):
        """
        calculates linear squared  distance between points, abs value of difference

        :param a: point A
        :type a: int
        :param b: point B
        :type b: int
        :return: distance
        :rtype: float
        """

        if b:
            tmp = a - b
        else:
            o = self.origin
            tmp = a - self.origin

        return np.abs(tmp)



if __name__ == "__main__":


    force_field = ForceField(linear=False)
    dst = force_field.get_eu_distance((1, 1, 1), (3, 3, 3))
    print(force_field.get_value()(np.sqrt(dst)))





