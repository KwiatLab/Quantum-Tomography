import unittest
import QuantumTomography as qLib
import numpy as np
import numpy.testing as tests

class MyTestCase(unittest.TestCase):

    # This tests phaserToComplex,complexToPhaser,removeGlobalPhase and stateToString
    def test_stateToString(self):
        phase1 = np.array([0.786236, 3.366])
        complex1 = qLib.phaserToComplex(phase1)
        phase1_2 = qLib.complexToPhaser(complex1)
        complex1_2 = qLib.phaserToComplex(phase1_2)
        tests.assert_array_almost_equal(phase1,phase1_2)
        tests.assert_almost_equal(complex1,complex1_2)

        complex2 = .5 + 5j
        phase2 = qLib.complexToPhaser(complex2)
        complex2_2 = qLib.phaserToComplex(phase2)
        phase2_2 = qLib.complexToPhaser(complex2_2)
        tests.assert_array_almost_equal(phase2,phase2_2)
        tests.assert_almost_equal(complex2,complex2_2)

        state = np.array([-.5 + .5j, -.5 - .5j])

        self.assertEqual('0.707|H> + i0.707|V>',qLib.stateToString(state))


if __name__ == '__main__':
    unittest.main()
