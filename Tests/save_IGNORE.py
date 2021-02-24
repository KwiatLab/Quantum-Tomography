from __future__ import print_function
from SaveRun import SaveRun
import unittest
import QuantumTomography as qLib
import warnings
warnings.filterwarnings("ignore")


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is UNTRACKED by git, and can be used however. 
By default it returns true. Most likely if you are running tests you'll want to save the
 results and use the SaveRun class, an example of using that class is in the save_RandomStates"""

class TestSum(unittest.TestCase):
    def test_main(self):
        t = qLib.Tomography()
        t.setConfSetting("NQubits",3)
        dataMatrix = t.getTomoInputTemplate()
        print(dataMatrix.shape)
        self.assertEqual(1, 1)

if __name__ == '__main__':
    unittest.main()
