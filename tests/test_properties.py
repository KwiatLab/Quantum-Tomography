import unittest
import quantum_tomography as qLib
import numpy as np
import numpy.testing as tests


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"Attention! These tests run on the version that your environment uses. See readme for details"

class Test_Properties(unittest.TestCase):
    # Test linear entropy and purity
    def test_purity(self):
        # Setup
        numStates = 1000
        for qDim in [1,2,3,4]:
            pure_states = [qLib.toDensity(qLib.random_pure_state(qDim)) for x in range(numStates)]
            puritys = np.array([qLib.purity(x) for x in pure_states])
            linear_entropys = np.array([qLib.linear_entropy(x) for x in pure_states])

            # Tests
            self.assertEqual(any(puritys<.95), False) # Purities of random pure states should be 1
            tests.assert_array_almost_equal(linear_entropys, 1-puritys) # Linear entropy = 1-purity

    # Test concurrence, tangle, and negativity
    def test_concurrence(self):
        # Setup
        bell_states = [qLib.toDensity(qLib.random_bell_state(2)) for x in range(100)]
        cons = np.array([qLib.concurrence(x) for x in bell_states])
        tangles = np.array([qLib.tangle(x) for x in bell_states])
        negs = np.array([qLib.negativity(x) for x in bell_states])
        # Tests
        self.assertEqual(any(cons < .95), False) # Concurrence of random bell states should be 1
        self.assertEqual(any(negs < .95), False)  # Negativity of random bell states should be 1
        tests.assert_array_almost_equal(tangles, cons**2) # tangle = concurrence**2

    # Test entropy
    def test_entropy(self):
        for qDim in [1,2,3,4]:
            # Setup
            bell_states = [qLib.toDensity(qLib.random_pure_state(qDim)) for x in range(1000)]
            entrops = np.array([qLib.entropy(x) for x in bell_states])
            # Tests
            self.assertEqual(any(entrops > .05), False)  # Entropy of pure states should be 0

    # Test fidelity
    def test_fidelity(self):
        for qDim in [1,2,3,4]:
            # Setup
            numStates = 1000
            pure_states = [qLib.random_pure_state(qDim) for x in range(numStates)]
            density_states = [qLib.random_density_state(qDim) for x in range(numStates)]
            bell_states = [qLib.toDensity(qLib.random_pure_state(qDim)) for x in range(numStates)]

            pure_states_almost = [qLib.performOperation(x,.9999999*np.eye(2**qDim)) for x in pure_states]
            density_states_almost = [qLib.performOperation(x,.9999999*np.eye(2**qDim)) for x in density_states]
            bell_states_almost = [qLib.performOperation(x,.9999999*np.eye(2**qDim)) for x in bell_states]

            fidels_pure_pp = []
            fidels_pure_dp = []
            fidels_pure_pd = []
            fidels_pure_dd = []
            fidels_density = []
            fidels_bell = []
            for i in range(numStates):
                fidels_pure_pp.append(qLib.fidelity(pure_states_almost[i],
                                                    pure_states[i]))
                fidels_pure_pd.append(qLib.fidelity(pure_states_almost[i],
                                                    qLib.toDensity(pure_states[i])))
                fidels_pure_dp.append(qLib.fidelity(qLib.toDensity(pure_states_almost[i]),
                                                    pure_states[i]))
                fidels_pure_dd.append(qLib.fidelity(qLib.toDensity(pure_states_almost[i]),
                                                    qLib.toDensity(pure_states[i])))
                fidels_density.append(qLib.fidelity(density_states_almost[i],
                                                    density_states[i]))
                fidels_bell.append(qLib.fidelity(bell_states_almost[i],
                                                    bell_states[i]))

            self.assertEqual(any(np.array(fidels_pure_pp) < .95), False)
            self.assertEqual(any(np.array(fidels_pure_pd) < .95), False)
            self.assertEqual(any(np.array(fidels_pure_dp) < .95), False)
            self.assertEqual(any(np.array(fidels_pure_dd) < .95), False)
            self.assertEqual(any(np.array(fidels_density) < .95), False)
            self.assertEqual(any(np.array(fidels_bell) < .95), False)
            randomFidels = np.array([qLib.fidelity(qLib.random_density_state(qDim),
                                                   qLib.random_density_state(qDim)) for x in range(numStates)])
            for i in range(int(2**qDim)):
                state_i = np.zeros(2**qDim)
                state_i[i] = 1
                for j in range(int(2 ** qDim)):
                    state_j = np.zeros(2**qDim)
                    state_j[j] = 1
                    if i == j:
                        self.assertAlmostEqual(qLib.fidelity(state_i,state_j),1)
                    else:
                        self.assertAlmostEqual(qLib.fidelity(state_i,state_j),0)



if __name__ == '__main__':
    unittest.main()
