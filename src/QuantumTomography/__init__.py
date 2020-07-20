
import sys
from .TomoClass import *
from .TomoDisplay import *
"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""



# # # # # # # # # #  #
# # TO THE DEVELOPER #
# # # # # # # # # #  #

# Everything in this code is based on the Matlab files you can download on our tomography website, a more detailed
# comments and definations of variables could be found in those m files, The function names are as same as possible.
#
# 'UseDerivatives' is used automatically when doing 1 detctor per qubit, and is not used when doing 2 detectors per qubit.
# This cause the program to run faster and provide more accurate results in both cases
