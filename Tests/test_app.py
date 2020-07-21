from __future__ import print_function
from TestRun import runTest
import QuantumTomography as qLib
import numpy as np
import unittest
import warnings
warnings.filterwarnings("ignore")


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

class TestSum(unittest.TestCase):
    #    1 Qubit
    def test_N1_e0_a0_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([1, 0, 0, 0, 0, 0, 0, 20]), 1 )
    def test_N1_e3_a0_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([1, 3, 0, 0, 0, 0, 0, 20]), 1 )
    def test_N1_e0_a0_d1_c0_b0_dr0(self):
        self.assertEqual(runTest([1, 0, 0, 1, 0, 0, 0, 20]), 1 )
    def test_N1_e0_a0_d0_c1_b0_dr0(self):
        self.assertEqual(runTest([1, 0, 0, 0, 1, 0, 0, 20]), 1 )
    def test_N1_e0_a0_d0_c0_b0_dr1(self):
        self.assertEqual(runTest([1, 0, 0, 0, 0, 1, 0, 20]), 1 )

    #     2 qubits
    def test_N2_e0_a0_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([2, 0, 0, 0, 0, 0, 0, 4]), 1 )
    def test_N2_e3_a0_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([2, 3, 0, 0, 0, 0, 0, 1]), 1 )
    def test_N2_e0_a1_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([2, 0, 1, 0, 0, 0, 0, 4]), 1 )
    def test_N2_e0_a0_d1_c0_b0_dr0(self):
        self.assertEqual(runTest([2, 0, 0, 1, 0, 0, 0, 4]), 1 )
    def test_N2_e0_a0_d0_c1_b0_dr0(self):
        self.assertEqual(runTest([2, 0, 0, 0, 1, 0, 0, 4]), 1 )
    def test_N2_e0_a0_d0_c0_b1_dr0(self):
        self.assertEqual(runTest([2, 0, 0, 0, 0, 1, 0, 4]), 1 )
    def test_N2_e0_a0_d0_c0_b0_dr1(self):
        self.assertEqual(runTest([2, 0, 0, 0, 0, 1, 0, 4]), 1 )
    def test_N3_e0_a0_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([3, 0, 0, 0, 0, 1, 0, 1]), 1 )
    def test_fidelity(self):
        for x in range(5):
            numQubits = np.random.randint(1, 4)
            state = np.random.beta(.5, .5, 2 ** numQubits) + 1j * np.random.beta(.5, .5, 2 ** numQubits)
            state = state / np.sqrt(np.dot(state, state.conj()))
            density = qLib.toDensity(state)
            self.assertGreater(qLib.fidelity(state,state), .95)
            self.assertGreater(qLib.fidelity(density,density), .95)
    def test_purity_linEntropy_pure(self):
        pureStates = np.zeros((6, 2), dtype=complex)
        pureStates[0] = np.array([1, 0], dtype=complex)
        pureStates[1] = np.array([0, 1], dtype=complex)
        pureStates[2] = np.array([.7071, .7071], dtype=complex)
        pureStates[3] = np.array([.7071, -.7071], dtype=complex)
        pureStates[4] = np.array([.7071, .7071j], dtype=complex)
        pureStates[5] = np.array([.7071, -.7071j], dtype=complex)
        for x in pureStates:
            state = qLib.toDensity(x)
            self.assertGreater(qLib.purity(state), .95)
            self.assertAlmostEquals(qLib.purity(state)+qLib.linear_entropy(state), 1)

    def test_purity_linEntropy_mixed(self):
        mixedStates = list()
        mixedStates.append(np.array([0, 0]))
        mixedStates.append(np.array([0, 0, 0, 0]))
        mixedStates.append(np.array([0, 0, 0, 0, 0, 0, 0, 0]))
        for x in mixedStates:
            state = qLib.toDensity(x)
            self.assertLess(qLib.purity(state), .05)
            self.assertAlmostEquals(qLib.purity(state) + qLib.linear_entropy(state), 1)

    def test_concurrence_tangle(self):
        bellStates = np.zeros((4, 4), dtype=complex)
        bellStates[0] = np.array([.7071, 0, 0, .7071], dtype=complex)
        bellStates[1] = np.array([.7071, 0, 0, -.7071], dtype=complex)
        bellStates[2] = np.array([0, .7071, .7071, 0], dtype=complex)
        bellStates[3] = np.array([0, .7071, -.7071, 0], dtype=complex)
        for x in bellStates:
            state = qLib.toDensity(x)
            self.assertGreater(qLib.concurrence(state), .95)
            self.assertGreater(qLib.tangle(state), .95)
import numpy as np
if __name__ == '__main__':
    unittest.main()