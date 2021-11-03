"""
tests for ForceField class
"""


from  sim.sim.force_field import  ForceField



def test_get_value():
    f_f = ForceField(amplitude=2)
    assert f_f.get_value()(1) == -2



def test_get_distance():
    f_f = ForceField()
    distance = f_f.get_eu_distance((1, 1, 1), (3, 3, 3))
    assert  distance == 12


def test_get_distance1():
    f_f = ForceField()
    distance = f_f.get_eu_distance((1, 1, 1))
    assert  distance == 3


def test_get_distance2():

    f_f = ForceField(linear=True)

    distance = f_f.get_distance(10)
    assert  distance == 10




def test_get_distance3():

    f_f = ForceField(linear=False)

    distance = f_f.get_distance((10,10,10))
    assert  distance == 300


def test_get_distance4():

    f_f = ForceField(linear=False)

    distance = f_f.get_distance((10,10,10), (20,20,20))
    assert  distance == 300


def test_get_distance5():
    f_f = ForceField(linear=True)

    distance = f_f.get_distance(40, 10)
    assert distance == 30











