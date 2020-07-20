from __future__ import print_function
import numpy as np

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""



# These are various helper functions used in other main functions in the TomoFunction file


# Helper function for partial Transpose in TomoFunctions
def partial_transpose_helper(m, d):
    if m.shape[0] == d:
        val = m.transpose()
    else:
        na = np.int(d)
        nb = np.int(len(m)/d)
        y = np.zeros([nb, nb, na, na])+0j
        val = np.zeros([len(m), len(m)])+0j
        for j in range(nb):
            for k in range(nb):
                y[j, k] = m[j*nb:j*nb+na, k*nb:k*nb+na]
        for j in range(nb):
            for k in range(nb):
                val[j*nb:j*nb+na, k*nb:k*nb+na] = y[k, j]

    return val

# performs the operation on the density matrix
def densityOperation(psi, gate):
    return np.matmul(np.matmul(gate, psi), np.conjugate(np.transpose(gate)))

# performs the operation on the ket state
def ketOperation(psi, gate):
    return np.matmul(gate, psi)
