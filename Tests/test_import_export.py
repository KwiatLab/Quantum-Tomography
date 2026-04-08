import unittest
from pathlib import Path
import QuantumTomography as qLib
from QuantumTomography.Utilities import OLD_FORMAT_CONFIG_KEYS
import numpy as np
import numpy.testing as tests
from TestRun import runTests

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"Attention! These tests run on the version that your environment uses. See readme for details"

TESTS_DIR = Path(__file__).parent
EXAMPLES_DIR = TESTS_DIR.parent / "ExampleFiles"

# Run tomographies
[[Tomo_Object_1, Fidelity_with_Original, Original_Purity, Total_Time]] = runTests(1, 1)
[[Tomo_Object_2, Fidelity_with_Original, Original_Purity, Total_Time]] = runTests(2, 1)
[[Tomo_Object_3, Fidelity_with_Original, Original_Purity, Total_Time]] = runTests(
    2,
    1,
    test2Det=True,
    testAccCorr=True,
    testCrossTalk=True,
    testDrift=True,
    testBell=True,
    errBounds=3,
    method="linear",
)


class Test_Import_Export(unittest.TestCase):
    def test_eval(self):
        filename = str(TESTS_DIR / "Test_States" / "rand_eval.txt")
        for Tomo_Object in [Tomo_Object_1, Tomo_Object_2, Tomo_Object_3]:
            # Export
            Tomo_Object.exportToEval(filename)

            # Import
            Tomo_Object_copy = qLib.Tomography()
            Tomo_Object_copy.importEval(filename)

            # Make sure exported conf settings are the same (only compare keys that are roundtripped)
            for k in OLD_FORMAT_CONFIG_KEYS:
                if not isinstance(Tomo_Object_copy.conf[k], np.ndarray):
                    self.assertEqual(Tomo_Object.conf[k], Tomo_Object_copy.conf[k])

            # make sure inputs are the same
            tests.assert_array_equal(Tomo_Object.last_input, Tomo_Object_copy.last_input)
            tests.assert_array_equal(Tomo_Object.intensities, Tomo_Object_copy.intensities)

            # Run tomographies and make sure estimates are the same
            tests.assert_array_equal(Tomo_Object.last_rho, Tomo_Object_copy.last_rho)
            self.assertEqual(Tomo_Object.last_intensity, Tomo_Object_copy.last_intensity)
            self.assertEqual(Tomo_Object.last_fval, Tomo_Object_copy.last_fval)

    def test_conf_and_data(self):
        filename_c = str(TESTS_DIR / "Test_States" / "rand_conf.txt")
        filename_d = str(TESTS_DIR / "Test_States" / "rand_data.txt")
        for Tomo_Object in [Tomo_Object_1, Tomo_Object_2, Tomo_Object_3]:
            # Export
            Tomo_Object.exportToConf(filename_c)
            Tomo_Object.exportToData(filename_d)

            # Import
            Tomo_Object_copy = qLib.Tomography()
            Tomo_Object_copy.importConf(filename_c)
            Tomo_Object_copy.importData(filename_d)

            # Make sure exported conf settings are the same (only compare keys that are roundtripped)
            for k in OLD_FORMAT_CONFIG_KEYS:
                if not isinstance(Tomo_Object_copy.conf[k], np.ndarray):
                    self.assertEqual(Tomo_Object.conf[k], Tomo_Object_copy.conf[k])

            # make sure inputs are the same
            tests.assert_array_equal(Tomo_Object.last_input, Tomo_Object_copy.last_input)
            tests.assert_array_equal(Tomo_Object.intensities, Tomo_Object_copy.intensities)

            # Run tomographies and make sure estimates are the same
            tests.assert_array_equal(Tomo_Object.last_rho, Tomo_Object_copy.last_rho)
            self.assertEqual(Tomo_Object.last_intensity, Tomo_Object_copy.last_intensity)
            self.assertEqual(Tomo_Object.last_fval, Tomo_Object_copy.last_fval)

    def test_printLastOutput(self):
        data_files = [
            EXAMPLES_DIR / "1_qubit_example.json",
            EXAMPLES_DIR / "bell_state_example.json",
            EXAMPLES_DIR / "2n_detector_example.json",
        ]
        for data_file in data_files:
            tomo = qLib.Tomography()
            tomo.import_data(str(data_file))
            tomo.run_tomography()
            tomo.printLastOutput()
            qLib.printLastOutput(tomo)

    def test_video_example(self):
        q = qLib.Tomography()
        q.import_data(str(EXAMPLES_DIR / "bell_state_example.json"))
        [rho_approx, intensity, fval] = q.run_tomography()


if __name__ == "__main__":
    unittest.main()
