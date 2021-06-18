import unittest
from TestRun import saveRunsFull,runTests

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

<<<<<<< HEAD
=======
class TestSum(unittest.TestCase):
    def test_main(self):
        self.assertEqual(runTest(1,10), 1)
        self.assertEqual(runTest(2,10), 1)
        self.assertEqual(runTest(1,1,method="BME"),1)
>>>>>>> 2f794aad6f0277c33d6271a52c2a30fdef8f3363

if __name__ == '__main__':
    # unittest.main()
    output = runTests(1, 2, method='MLE')