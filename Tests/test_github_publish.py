from __future__ import print_function
from gen_Full import runFull
from test_properties import Test_Properties
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

"""This test script is the one used by github when a new version is created. This will automatically
run and the results can be see in the actions tab"""

"Attention! These tests run on the version that your environment uses. see readme for details"

class TestFull(unittest.TestCase):
    #    1 Qubit
    def test_FULL(self):
        numErrors = runFull(nStates=1,saveStates=False)
        self.assertEqual(numErrors,0)



if __name__ == '__main__':
    unittest.TestLoader().loadTestsFromTestCase(Test_Properties)
    unittest.main()
