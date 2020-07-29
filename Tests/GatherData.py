from SaveTestRun import TestRun
import numpy as np

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is used to run tomography and save the results to the Results Folder.
This uses the SaveTestRun.py script which saves the data."""


for x in range(1,2):
    q = TestRun([x, 0, False, False, False, False, False, False])
    q.run()
