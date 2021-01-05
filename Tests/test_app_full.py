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

"""This script is intended to be a full test that is larger than the test_github.py. 
Currently does not run tests and returns true.
test_github_*.py should be the main testing script. Until that becomes too big (or some other reason that you want some tests not on github)
this script should remain unedited."""

class TestSum(unittest.TestCase):
    def test_main(self):
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()
