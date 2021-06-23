import unittest
import QuantumTomography as qLib
import numpy as np
import numpy.testing as tests
from TestRun import runTests

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"Attention! These tests run on the version that your environment uses. See readme for details"

# todo : add a test for cross talk. When doing tomography with a 50:50 BS you should get a compltely mixed state

# add a test where the error bounds is not specified, but then getProperties(5) is run after
# then do getProperties(3) is run and should not increase the number of monte carlo states.
# assert that the number of monte carlo states is 5 after both getProperties

class Test_Properties(unittest.TestCase):

    # # Error bounds is not specified, then getProperties(5) is run after
    # # then do getProperties(3) is run and should not increase the number of monte carlo states.
    # def test_getProperties(self):
    #     for i in [1,2]:
    #         tomo = runTests(numQubits=i, nStates=1)[0][0]
    #         self.assertEqual(len(tomo.mont_carlo_states),1)
    #         tomo.getProperties(5)
    #         self.assertEqual(len(tomo.mont_carlo_states), 6)
    #         tomo.getProperties(3)
    #         self.assertEqual(len(tomo.mont_carlo_states), 6)

    def test_getProperties(self):
        cross = np.ones(1,dtype=complex)
        for i in [1,2,3]:
            # When doing tomography with a BS that gets rid of all information you should
            # get a completely mixed state
            cross = np.kron(1/2*np.array([[1,1],[1,1]]),cross)
            list_of_tomo = runTests(numQubits=i, nStates=10,testCrossTalk=cross)
            for t in list_of_tomo:
                print(t[0].last_rho)
                # tests.assert_array_almost_equal(t[0].last_rho,1/2**i*np.ones((2**i,2**i)))

if __name__ == '__main__':
    unittest.main()
