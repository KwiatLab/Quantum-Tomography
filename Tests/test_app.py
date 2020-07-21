from __future__ import print_function
from TestRun import runTest


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""


def simple_general_test():

    assert(runTest([1, 0, 0, 0, 0, 0, 0, 50]) == 1 )
    assert(runTest([1, 3, 0, 0, 0, 0, 0, 50]) == 1 )
    assert(runTest([1, 0, 0, 1, 0, 0, 0, 50]) == 1 )
    # assert(runTest([1, 0, 0, 0, 1, 0, 0, 50]) == 1 )
    # assert(runTest([1, 0, 0, 0, 0, 0, 1, 50]) == 1 )
    #
    # assert(runTest([2, 0, 0, 0, 0, 0, 0, 50]) == 1 )
    # assert(runTest([2, 3, 0, 0, 0, 0, 0, 50]) == 1 )
    # assert(runTest([2, 0, 1, 0, 0, 0, 0, 50]) == 1 )
    # assert(runTest([2, 0, 0, 1, 0, 0, 0, 50]) == 1 )
    # assert(runTest([2, 0, 0, 0, 1, 0, 0, 50]) == 1 )
    # assert(runTest([2, 0, 0, 0, 0, 1, 0, 50]) == 1 )
    # assert(runTest([2, 0, 0, 0, 0, 0, 1, 50]) == 1 )
    #
    # assert(runTest([3, 0, 0, 0, 0, 0, 0, 50]) == 1 )
    
simple_general_test()
