class  CubicCell:
    """
    represents confining cell for the polymer
    """

    def __init__(self, a,b,c):
        """
        init the cell

        :param a: ``OX`` dimension
        :type a: int
        :param b:  ``OY`` dimension
        :type b: int
        :param c: ``OZ`` dimension
        :type c:  int
        """

        self.A = a
        self.B = b
        self.C = c
        self.center = self.set_center()


    def __str__(self):
        return "the cell is a parallelogram with the size (X,Y,Z) = (%i, %i, %i)"%(self.A, self.B, self.C)


    def set_center(self):
        """
        calculates the center of the cell

        :return: list of 3 numbers
        :rtype: list
        """
        return [self.A//2, self.B//2, self.C//2]


    def get_center(self):
        """
        returns the center of the  cell

        :return: list of integers
        :rtype: list
        """
        return  self.center
