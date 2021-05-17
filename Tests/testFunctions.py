import unittest
import QuantumTomography as qlib
import numpy as np
import numpy.linalg as la


class functionTesting(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(functionTesting, self).__init__(*args, **kwargs)
        self.tomoinp1 = np.array([[1,0,1,1,0],[1,0,1,0,1],[1,0,1,0.7071,0.7071],[1,0,1,0.7071,-0.7071],[1,0,1,0.7071,0.7071j],[1,0,1,0.7071,-0.7071j]])
        self.tomoinp2 = np.array([[1, 0, 0, 500, 1, 0, 1, 0], [1, 0, 0, 500, 1, 0, 0, 1], [1, 0, 0, 500, 1, 0, 0.7071, 0.7071],
                             [1, 0, 0, 500, 1, 0, 0.7071, 0.7071j], [1, 0, 0, 500, 0, 1, 1, 0], [1, 0, 0, 500, 0, 1, 0, 1],
                             [1, 0, 0, 500, 0, 1, 0.7071, 0.7071], [1, 0, 0, 500, 0, 1, 0.7071, 0.7071j],
                             [1, 0, 0, 500, 0.7071, 0.7071, 1, 0], [1, 0, 0, 500, 0.7071, 0.7071, 0, 1],
                             [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071], [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071j],
                             [1, 0, 0, 500, 0.7071, 0.7071j, 1, 0], [1, 0, 0, 500, 0.7071, 0.7071j, 0, 1],
                             [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071], [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071j]])
        self.tomoObj1 = qlib.Tomography(1)
        self.tomoObj2 = qlib.Tomography()

    def test_getCoincidences(self):
        #setting up new inputs to be changed
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        #testing different random measurement values
        for i in range(500, 1500, 200):
            rand_coincidence1 = np.random.binomial(i, 0.5, 6)
            rand_coincidence2 = np.random.binomial(i, 0.5, 16)
            curinp1[:,2] = rand_coincidence1
            curinp2[:,3] = rand_coincidence2
            tomo = self.tomoObj1.state_tomography(curinp1)
            tomo2 = self.tomoObj2.state_tomography(curinp2)
            self.assertTrue((self.tomoObj1.getCoincidences() == rand_coincidence1).all())
            self.assertTrue((self.tomoObj2.getCoincidences() == rand_coincidence2).all())


    def test_getSingles(self):
        #setting up new inputs to be changed
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        #test case with 0 everywhere
        tomo = self.tomoObj1.state_tomography(curinp1)
        tomo2 = self.tomoObj2.state_tomography(curinp2)
        self.assertTrue((self.tomoObj1.getSingles() == np.zeros(6)).all())
        self.assertTrue((self.tomoObj2.getSingles() == np.zeros((16,2))).all())

        #test cases with different random singles values
        for i in range(1, 20, 4):
            singles_onebit = np.random.binomial(i, 0.5, 6)
            singles_twobit = np.random.binomial(i, 0.5, (16,2))
            curinp1[:, 1] = singles_onebit
            curinp2[:, 1:3] = singles_twobit
            tomo = self.tomoObj1.state_tomography(curinp1)
            tomo2 = self.tomoObj2.state_tomography(curinp2)
            self.assertTrue((np.real(np.transpose(self.tomoObj1.getSingles()))[0] == singles_onebit).all())
            self.assertTrue((np.real(self.tomoObj2.getSingles()) == singles_twobit).all())


    def test_getTimes(self):
        #make copies of inputs to change
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        #test case with all ones for time
        tomo  = self.tomoObj1.state_tomography(curinp1)
        tomo2 = self.tomoObj2.state_tomography(curinp2)
        self.assertTrue((self.tomoObj1.getTimes() == np.ones(6)).all())
        self.assertTrue((self.tomoObj2.getTimes() == np.ones(16)).all())

        #testing some different cases:
        for i in range(3):
            randTime1 = np.random.randint(1, 10, 6)
            randTime2 = np.random.randint(1, 10, 16)
            curinp1[:, 0] = randTime1
            curinp2[:, 0] = randTime2

            tomo  = self.tomoObj1.state_tomography(curinp1)
            tomo2 = self.tomoObj2.state_tomography(curinp2)
            self.assertEqual((self.tomoObj1.getTimes() == randTime1).all())
            self.assertEqual((self.tomoObj2.getTime() == randTime2).all())


if __name__ == '__main__':
    unittest.main()
