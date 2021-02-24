import numpy as np
import sys
sys.path.append("../Quantum-Tomography_Git")
import QuantumTomography as qLib
from SaveRun import SaveRun

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is used to run tomography against the a specific set of states. The measurements and couts 
(aka the data) used for the tomographies are in the Test_States Foulder. standard_test_data, and save the Results"""

"Attention! These tests run on the version that your environment uses. see readme for details"

test1 = SaveRun([1, 0, 0, 1, 0, 0, 0, 25])
test1.run()