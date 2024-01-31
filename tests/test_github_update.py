from __future__ import print_function
from TestRun import runTests
from test_properties import Test_Properties
from test_functions import Test_Functions
import unittest
import warnings
warnings.filterwarnings("ignore")


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"""This test script is the one used by github when there is a pull request or push. This will automatically
run and the results can be see in the actions tab"""

"Attention! These tests run on the version that your environment uses. see readme for details"

class TestQuick(unittest.TestCase):
    #    1 Qubit
    def test_N1_e0_a0_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(1,50)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N1_e3_a0_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(1,20,errBounds=3)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N1_e0_a0_d1_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(1,20,test2Det=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N1_e0_a0_d0_c1_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(1,20,testCrossTalk=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N1_e0_a0_d0_c0_b0_dr1(self):
        numTotalErrors = 0
        Tomographys = runTests(1,20,testDrift=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)

    # 2 qubits
    def test_N2_e0_a0_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(2,15)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N2_e3_a0_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(2,5,errBounds=3)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N2_e0_a1_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(2,5,testAccCorr=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N2_e0_a0_d1_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(2,5,test2Det=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N2_e0_a0_d0_c1_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(2,5,testCrossTalk=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N2_e0_a0_d0_c0_b1_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(2,5,testBell=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)
    def test_N2_e0_a0_d0_c0_b0_dr1(self):
        numTotalErrors = 0
        Tomographys = runTests(2,5,testDrift=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors,0)

    # 3 qubits
    def test_N3_e0_a0_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(3, 1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors, 0)

    def test_N3_e3_a0_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(3, 1, errBounds=3)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors, 0)

    def test_N3_e0_a1_d0_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(3, 1, testAccCorr=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors, 0)

    def test_N3_e0_a0_d1_c0_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(3, 1, test2Det=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors, 0)

    def test_N3_e0_a0_d0_c1_b0_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(3, 1, testCrossTalk=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors, 0)

    def test_N3_e0_a0_d0_c0_b1_dr0(self):
        numTotalErrors = 0
        Tomographys = runTests(3, 1, testBell=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors, 0)

    def test_N3_e0_a0_d0_c0_b0_dr1(self):
        numTotalErrors = 0
        Tomographys = runTests(3, 1, testDrift=1)
        for t in Tomographys:
            [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
            if Original_Purity < 0:
                numTotalErrors += 1
        self.assertEqual(numTotalErrors, 0)
    

if __name__ == '__main__':
    unittest.TestLoader().loadTestsFromTestCase(Test_Properties)
    unittest.TestLoader().loadTestsFromTestCase(Test_Functions)
    unittest.main()