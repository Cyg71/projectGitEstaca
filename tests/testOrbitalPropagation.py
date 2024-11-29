from orbitalPropagation import *

def test_Rot():
    assert np.linalg.det(Rot(0, 0)) == 1
    assert np.linalg.det(Rot(0, 1)) == 1
    assert np.linalg.det(Rot(0, 2)) == 1

def test_radius():
    a = 6900e3
    assert radius(a, 0, 0) == a