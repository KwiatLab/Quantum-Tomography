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

"ATTENTION! This script may or may not be testing the src code. see readme"

test1 = SaveRun([1, 0, 0, 0, 0, 0, 0, 20])
test1.run()