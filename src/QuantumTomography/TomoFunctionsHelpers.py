from __future__ import print_function
import numpy as np
import numpy.random as rand

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


# These are various helper functions used in other main functions in the TomoFunction file


"""
    random_ginibre(D)
    Desc: Returns a random matrix from the Ginibre ensemble of size DxD. 
    This is a complex matrix whos elements are a+ib | a,b iid. Norm(0,1)

    Parameters
    ----------
    D :int
        The dimension of the Matrix

    Returns
    -------
    mat : ndarray with shape = (2^N, 2^N)
        The random matrix
    """


def random_ginibre(D=2):
    mat = np.zeros((D, D), dtype=complex)
    for i in range(D):
        for j in range(D):
            mat[i, j] = rand.normal(0, 1) + rand.normal(0, 1) * 1j

    return mat


def phaserToComplex(phaser):
    Magn = phaser[0]
    Phase = phaser[1]
    complex = Magn * (np.cos(Phase) + 1j * np.sin(Phase))
    return complex


def complexToPhaser(complex):
    Magn = np.absolute(complex)
    Phase = np.angle(complex)
    if Phase < 0:
        Phase += 2 * np.pi
    return np.array([Magn, Phase])


def isStateVector(state):
    return np.ndim(state) == 1
