import QuantumTomography as qLib
import numpy as np
import unittest

class Test_Properties(unittest.TestCase):
    def test_example_2q2d(self):
        # Coincidence counts for 9 measurements.
        # Each row is the coincidence counts for detector pairs 1,2,3 and 4 respectfully
        coincidence_counts = np.array([[497, 0, 0, 503],
                                       [256, 265, 250, 229],
                                       [235, 268, 242, 255],
                                       [254, 262, 249, 235],
                                       [521, 0, 0, 479],
                                       [235, 268, 248, 249],
                                       [265, 229, 280, 226],
                                       [253, 242, 247, 258],
                                       [0, 495, 505, 0]])

        # Measurement basis for 9 measurements
        # In each row the first two (possible complex) numbers alpha and beta represent the state that the first qubit
        # is projected onto when it ends up at detector 1.
        # The next two numbers is the state the second qubit is projected onto when it ends up at detector 2.
        measurements = np.array([[1 + 0j, 0j, 1 + 0j, 0j],
                                 [1 + 0j, 0j, 0.70710678 + 0j, 0.7071067811865476 + 0j],
                                 [1 + 0j, 0j, 0.70710678 + 0j, 0.70710678j],
                                 [0.70710678 + 0j, 0.70710678 + 0j, 1 + 0j, 0j],
                                 [0.70710678 + 0j, 0.70710678 + 0j, 0.70710678 + 0j, 0.70710678 + 0j],
                                 [0.70710678 + 0j, 0.70710678 + 0j, 0.70710678 + 0j, 0.70710678j],
                                 [0.70710678 + 0j, 0.70710678j, 1 + 0j, 0j],
                                 [0.70710678 + 0j, 0.70710678j, 0.70710678 + 0j, 0.70710678 + 0j],
                                 [0.70710678 + 0j, 0.70710678j, 0.70710678 + 0j, 0.7071067811865476j], ])

        # Initiate tomography object
        tomo_obj = qLib.Tomography()

        # Run tomography
        [rho_approx, intensity, fval] = tomo_obj.StateTomography(measurements, coincidence_counts)

        # Print Results
        tomo_obj.printLastOutput()
        print('----------------')
        bell_state = 1 / np.sqrt(2) * np.array([1, 0, 0, 1], dtype=complex)
        print('Fidelity with actual : ' + str(qLib.fidelity(bell_state, rho_approx)))
        self.assertTrue(1 - qLib.fidelity(bell_state, rho_approx) < 10-6)


if __name__ == '__main__':
    unittest.main()
