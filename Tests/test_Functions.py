import unittest
import QuantumTomography as qlib
import numpy as np
import warnings
warnings.filterwarnings("ignore")

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Shallat'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"Attention! These tests run on the version that your environment uses. see readme for details"

class FunctionTesting(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(FunctionTesting, self).__init__(*args, **kwargs)

        # initializing the inputs to be used in the following test tomographies
        self.tomoinp1 = np.array([[1,0,500,1,0],[1,0,500,0,1],[1,0,500,0.7071,0.7071],[1,0,500,0.7071,-0.7071],[1,0,500,0.7071,0.7071j],[1,0,500,0.7071,-0.7071j]])
        self.tomoinp2 = np.array([[1, 0, 0, 500, 1, 0, 1, 0], [1, 0, 0, 500, 1, 0, 0, 1], [1, 0, 0, 500, 1, 0, 0.7071, 0.7071],
                             [1, 0, 0, 500, 1, 0, 0.7071, 0.7071j], [1, 0, 0, 500, 0, 1, 1, 0], [1, 0, 0, 500, 0, 1, 0, 1],
                             [1, 0, 0, 500, 0, 1, 0.7071, 0.7071], [1, 0, 0, 500, 0, 1, 0.7071, 0.7071j],
                             [1, 0, 0, 500, 0.7071, 0.7071, 1, 0], [1, 0, 0, 500, 0.7071, 0.7071, 0, 1],
                             [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071], [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071j],
                             [1, 0, 0, 500, 0.7071, 0.7071j, 1, 0], [1, 0, 0, 500, 0.7071, 0.7071j, 0, 1],
                             [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071], [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071j]])


    def test_getCoincidences(self):
        # setting up tomography objects
        tomoObj1 = qlib.Tomography(1)
        tomoObj2 = qlib.Tomography(2)
        # setting up new inputs to be changed
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        # testing different random measurement values
        for i in range(500, 1500, 200):
            rand_coincidence1 = np.random.binomial(i, 0.5, 6)
            rand_coincidence2 = np.random.binomial(i, 0.5, 16)
            curinp1[:,2] = rand_coincidence1
            curinp2[:,3] = rand_coincidence2
            tomo = tomoObj1.state_tomography(curinp1)
            tomo2 = tomoObj2.state_tomography(curinp2)
            self.assertTrue((tomoObj1.getCoincidences() == rand_coincidence1).all())
            self.assertTrue((tomoObj2.getCoincidences() == rand_coincidence2).all())


    def test_getSingles(self):
        # setting up tomography objects
        tomoObj1 = qlib.Tomography(1)
        tomoObj2 = qlib.Tomography(2)
        # setting up new inputs to be changed
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        # test case with 0 everywhere
        tomo = tomoObj1.state_tomography(curinp1)
        tomo2 = tomoObj2.state_tomography(curinp2)
        self.assertTrue((tomoObj1.getSingles() == np.zeros(6)).all())
        self.assertTrue((tomoObj2.getSingles() == np.zeros((16,2))).all())

        # test cases with different random singles values
        for i in range(1, 20, 4):
            singles_onebit = np.random.binomial(i, 0.5, 6)
            singles_twobit = np.random.binomial(i, 0.5, (16,2))
            curinp1[:, 1] = singles_onebit
            curinp2[:, 1:3] = singles_twobit
            tomo = tomoObj1.state_tomography(curinp1)
            tomo2 = tomoObj2.state_tomography(curinp2)
            self.assertTrue((np.real(np.transpose(tomoObj1.getSingles()))[0] == singles_onebit).all())
            self.assertTrue((np.real(tomoObj2.getSingles()) == singles_twobit).all())


    def test_getTimes(self):
        # setting up tomography objects
        tomoObj1 = qlib.Tomography(1)
        tomoObj2 = qlib.Tomography(2)
        # make copies of inputs to change
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        # test case with all ones for time
        tomo  = tomoObj1.state_tomography(curinp1)
        tomo2 = tomoObj2.state_tomography(curinp2)
        self.assertTrue((tomoObj1.getTimes() == np.ones(6)).all())
        self.assertTrue((tomoObj2.getTimes() == np.ones(16)).all())

        # testing some different cases:
        for i in range(3):
            randTime1 = np.random.randint(1, 10, 6)
            randTime2 = np.random.randint(1, 10, 16)
            curinp1[:, 0] = randTime1
            curinp2[:, 0] = randTime2

            tomo  = tomoObj1.state_tomography(curinp1)
            tomo2 = tomoObj2.state_tomography(curinp2)
            self.assertTrue((tomoObj1.getTimes() == randTime1).all())
            self.assertTrue((tomoObj2.getTimes() == randTime2).all())

    def test_buildTomoInput(self):
        tomoObj1 = qlib.Tomography(1)
        tomoObj2 = qlib.Tomography(2)

        [times1, singles1, coincidences1, measurements11, measurements12] = self.tomoinp1.transpose()



        measurements1 = np.array([measurements11,measurements12]).transpose()
        built_input1 = tomoObj1.buildTomoInput(measurements1,coincidences1,crosstalk=-1,efficiency=0,time=times1,singles=singles1,window=0,error=0)
        self.assertTrue((built_input1==self.tomoinp1).all())


        [times2, singles21, singles22, coincidences2, measurements211, measurements212, measurements221, measurements222] = self.tomoinp2.transpose()
        singles2 = np.array([singles21, singles22]).transpose()

        measurements2 = np.array([measurements211,measurements212,measurements221,measurements222]).transpose()
        built_input2 = tomoObj2.buildTomoInput(measurements2, coincidences2, crosstalk=-1,efficiency=0,time=times2,singles=singles2,window=0,error=0)

        self.assertTrue((built_input2 == self.tomoinp2).all())

if __name__ == '__main__':
    unittest.main()
