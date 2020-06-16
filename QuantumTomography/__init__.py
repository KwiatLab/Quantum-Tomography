"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""
import sys

from .TomoClass import *
from .TomoDisplay import *

__author__ = 'Quoleon/Turro'
#################
## TO THE USER ##
#################

# To use the code to calculate tomography and to estimate error of a measurement, you can skip these tedious definations
# of functions and scroll down to the bottom of it, where marked as ##TESTING DATA##. There you can use sample data or
# import your own data. ##MAIN FUNCTION## is where the main method state_tomography() is called. This is also where the
# method displayOutput() is called which creates a graph of the density matrix and prints the error values
#What's new:
#
# 1. State Tomography can deal with n qubits now, by changing 'NQubits' option and input matrix. Error estimation and
#    Bell state searching are NOT working for N qubits.
#
# 2. Two types of random state generator is added, providing state with different purity.
#
# 3. Coincidence generator and tomo_input generator is added, to estimate how good a measurement can be.
#
# 4. Hedged MLE method is added, still testing.
#
# 5. Cross talk matrix generator is added, providing an easier way to input the property of beam splitter.
#
# 6. Faster MLE process.
#
#
# How to use:
#
# pythoneval.txt should contends following input data:
#
# tomo_input: array like, demension = 2.
#
# For n detectors:
# tomo_input[:, 0]: times
# tomo_input[:, np.arange(1, n_qubit+1)]: singles
# input[:, n_qubit+1]: coincidences
# input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
#
# For 2n detectors:
# tomo_input[:, 0]: times
# tomo_input[:, np.arange(1, 2*n_qubit+1)]: singles
# input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
# input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements
#
# intensity: array like, demension = 1, length = length of tomo_input
#
# conf['NQubits']: >=1, it will take a much longer time for more qubits.
# conf['NDetectors']: 1 or 2
# conf['ctalk']: [[C0->0, C0->1],[C1->0, C1->1]]
# conf['UseDerivative']: 0 or 1
# conf['Bellstate']: 0 or 1
# conf['DoErrorEstimation']: 0 or 1
# conf['DoDriftCorrection'] = 'no' or 'yes'
# conf['Window']: 0 or array like, demension = 1
# conf['Efficiency']: 0 or array like, demension = 1
# conf['Beta']: 0 to 0.5, depending on purity of state and total number of measurements.

"""Make sure to add numpy, scipy and matplotlib to your project interpreter before starting!"""

######################
## TO THE DEVELOPER ##
######################

# Everything in this code is based on the Matlab files you can download on our tomography website, a more detailed
# comments and definations of variables could be found in those m files, The function names are as same as possible.

# If you are debugging this program there is a section underneath the sample sets that involve some usefull debugging
# data sets. You can set the state or choose to have it be random and it will create data.
#
# 'UseDerivatives' is used automatically when doing 1 detctor per qubit, and is not used when doing 2 detectors per qubit.
# This cause the program to run faster and provide more accurate results in both cases
#
#
#   Drift Correction is just set off by default. and there is no implementation of it in this code.
#   I dont know the equation to properly implement it
