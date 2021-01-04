import numpy as np
import sys
sys.path.append("../Quantum-Tomography_Git")
import QuantumTomography as qLib
from SaveRun import TestRun

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is used to run tomo against the standard_test_data, and save the Results"""

"WARNING! These tests run on the published library installed in your pip version, not the code in the local directory."

test1 = TestRun([2, 0, 0, 0, 0, 0, 0, 20])
test1.run()