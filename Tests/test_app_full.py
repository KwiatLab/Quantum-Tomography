from __future__ import print_function
from TestRun import runTest
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
    def test_N1_e0_a0_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([1, 0, 0, 0, 0, 0, 0, 20]), 1)

    def test_N1_e3_a0_d0_c0_b0_dr0(self):
        self.assertEqual(runTest([1, 3, 0, 0, 0, 0, 0, 20]), 1)

if __name__ == '__main__':
    unittest.main()