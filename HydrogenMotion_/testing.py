import numpy as np
import unittest
import atom
import readfiles as rf
import motion
import subprocess as sp


class TEST():
    
    def __init__(self, vals):
        self.x=vals[0]
        self.y=vals[1]
        self.neighbors=None
    
    def set_neighbors(self, ar):
        self.neighbors = ar

#=====================================
class Test_atom(unittest.TestCase):
    
    def test_distance(self):
        a1=atom.ATOM([216, 1, 1.48, 13.2, 9.654])
        a2=atom.ATOM([134, 2, 16.1, 9.57, 9.567])
        atom.ATOM.set_box([[0.0, 16.29], [0.0, 16.29], [0.0, 16.29]])
        result = atom.distance(a1,a2)
        print(result)
    
    def test_dist_moved(self):
        a1=atom.ATOM([216, 1, 1.48, 13.2, 9.654])
        a2=atom.ATOM([134, 2, 16.1, 9.57, 9.567])
        atom.ATOM.set_box([[0.0, 16.29], [0.0, 16.29], [0.0, 16.29]])
        a1.coords=a2.coords
        result = a1.distance_moved()
        print(result)
#====================================
v1 = [1.5, 3.0]
v2 = [3.55, 2.22]
v3 = [2.7, 0.03]
t1=TEST(v1)
t2=TEST(v2)
t3=TEST(v3)
t1.set_neighbors(np.array([t2,t3]))
t2.set_neighbors(np.array([t1,t3]))
t3.set_neighbors(np.array([t2,t1]))
x=np.isin(t1,t2.neighbors)
x2="this is a test"
emp=[True, False, False, True, True]
mask=np.array([True, True, False, False, True])
x=100
sp.run(["mkdir", "/Users/diggs/Desktop/out_2/"])
if __name__ =='__main__':
    unittest.main()











