from __future__ import print_function
import scipy as sp
import numpy as np
from .TomoFunctionsHelpers import *

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""


# # # # # # # # # # # # #
# TOMOGRAPHY CALCULATE  #
# # # # # # # # # # # # #

#
# def i2array(i, ii, n):
#     nn = np.int(np.ceil((np.log(ii)/np.log(n))))
#     rv = np.zeros(nn)
#     for j in range(nn):
#         rv[j] = i/(n**(nn-j-1))
#         i % = n**(nn-j-1)
#     return rv
#
# # returns the tensor product of the two states
# def tensor_product(A, B):
#     a = np.ndim(A)
#     b = np.ndim(B)
#     if (a == 2) & (b == 2):
#         [n11, n12] = np.shape(A)
#         [n21, n22] = np.shape(B)
#         jj = n11 * n21
#         kk = n12 * n22
#         rv = np.zeros([jj, kk]) + 0j
#         for j in range(jj):
#             for k in range(kk):
#                 rv[j, k] = A[int(np.floor(j / n21))][int(np.floor(k / n22))] * B[j % n21][k % n22]
#     elif (a == 2) & (b == 1):
#         [n11, n12] = np.shape(A)
#         n21 = len(B)
#         jj = n11 * n21
#         kk = n12
#         rv = np.zeros([jj, kk]) + 0j
#         for j in range(jj):
#             for k in range(kk):
#                 rv[j, k] = A[int(np.floor(j / n21))][k] * B[j % n21]
#     elif (a == 1) & (b == 2):
#         [n21, n22] = np.shape(B)
#         n11 = len(A)
#         jj = n11 * n21
#         kk = n22
#         rv = np.zeros([jj, kk]) + 0j
#         for j in range(jj):
#             for k in range(kk):
#                 rv[j, k] = A[int(np.floor(j / n21))] * B[j % n21][ k]
#     elif (a == 1) & (b == 1):
#         n11 = len(A)
#         n21 = len(B)
#         jj = n11 * n21
#         rv = np.zeros(jj) + 0j
#         for j in range(jj):
#             rv[j] = A[int(np.floor(j / n21))] * B[j % n21]
#     elif (a == 0) | (b == 0):
#         rv = A * B
#
#     return rv
#
# def trace_dist(rho1, rho2):
#     # didn't checked, and would not be called in this version.
#     s1 = rho2stokes(rho1)
#     s2 = rho2stokes(rho2)
#     s = s1 - s2
#     val = np.sqrt(np.dot(s.conj().transpose(), s))/2
#
#     return val


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : ndarray
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def make_positive(rhog_in):
    d, v = np.linalg.eig(rhog_in)
    rhog = np.zeros(rhog_in.shape)
    for j in range(len(d)):
        rhog = rhog + np.abs(d[j])*np.outer(v[:, j], v[:, j].conj().transpose())
    rhog = (rhog + rhog.conj().transpose())/2.0

    return rhog


"""
    density2tm(rhog)
    Desc: Converts a density matrix into a lower t matrix.

    Parameters
    ----------
    rhog : ndarray
        Array to be converted.

    Returns
    -------
    tm : ndarray
        Lower t matrix defining the input matrix.
    """
def density2tm(rhog):
    d = rhog.shape[0]
    if d == 1:
        tm = np.real(np.sqrt(rhog))
        return tm

    tm = np.zeros(rhog.shape)+0j
    last_element = rhog[d-1][d-1]
    tm[d-1][d-1] = np.real(np.sqrt(last_element))
    if last_element > -.00000001:
        temp = rhog[d-1][0:(d-1)]
        tm[d-1][0:(d-1)] = temp/np.sqrt(last_element)
        # switched order of temp and temp.conj and transpose()
        recurse = np.hsplit(rhog[0:(d-1)], [d-1, d])[0] - np.outer(temp.conj().transpose(), temp)/last_element
    else:
        tm[d-1][0:(d-1)] = np.zeros(d)
        recurse = np.hsplit(rhog[0:(d-1)], [d-1, d])[0]
    for i in range(d-1):
        tm[i][0:(d-1)] = density2tm(recurse)[i][0:(d-1)]

    return tm

"""
    density2t(rhog)
    Desc: Converts a density matrix into a list of t values

    Parameters
    ----------
    rhog : ndarray
        Array to be converted.

    Returns
    -------
    t : ndarray
        List of t values defining the input matrix.
    """
def density2t(rhog):
    tm = density2tm(rhog)
    d = len(tm)

    idx = 0
    cur_length = d
    t = np.zeros(d**2)

    for j in range(d):
        t[np.arange(idx, idx+cur_length)] = np.real(np.diag(tm, -j))
        idx = idx + cur_length
        if j > 0:
            t[np.arange(idx, idx+cur_length)] = np.imag(np.diag(tm, -j))
            idx = idx + cur_length
        cur_length -= 1

    return t

"""
    toDensity(psiMat)
    Desc: Converts a pure state into a density matrix.

    Parameters
    ----------
    psiMat : ndarray
        Pure state to be converted.

    Returns
    -------
    rhog: ndarray
        Density Matrix of the input state.
    """
def toDensity(psiMat):
    return np.outer(psiMat, psiMat.conj())

#
# def one_in(idx, length):
#     val = np.zeros(length)
#     val[idx] = 1
#
#     return val


"""
    t_matrix(t)
    Desc: Converts a list of t values to an lower t matrix.

    Parameters
    ----------
    t : ndarray
        List of t values converted.

    Returns
    -------
    tm : ndarray
        Lower t matrix.
    """
def t_matrix(t):
    d = np.int(np.sqrt(len(t)))

    idx = 0
    cur_length = d
    tm = np.zeros([d, d])

    for j in range(np.int(d)):
        tm = tm + 1*np.diag(t[np.arange(idx, idx+cur_length)], -j)
        idx = idx + cur_length

        if j > 0:
            tm = tm + 1j*np.diag(t[np.arange(idx, idx+cur_length)], -j)
            idx = idx + cur_length

        cur_length -= 1

    return tm

"""
    t_to_density(t)
    Desc: Converts a list of t values to a density matrix.

    Parameters
    ----------
    t : ndarray
        List of t values converted.

    Returns
    -------
    rhog : ndarray
        Density Matrix.
    """
def t_to_density(t):
    tm = t_matrix(t)
    rhog = np.dot(tm.conj().transpose(), tm)

    return rhog


# # # # # # # # # #
# ERROR ESTIMATE  #
# # # # # # # # # #


"""
    fidelity(state1, state2)
    Desc: Calculates the fidelity between the two input states.

    Parameters
    ----------
    state1, state2 : ndarray
        Input arrays to calculate the fidelity between. Can be pure states or density matrices.

    Returns
    -------
    val : float
        The calculated fidelity.
    """
def fidelity(state1, state2):
    pure = 0
    rho1 = state1
    rho2 = state2
    if np.ndim(state1) == 1:
        state1 = toDensity(state1)
        rho1 = np.dot(state1, state1.conj().transpose())
        pure = 1
    elif state1.shape[1] == state1.shape[0]:
        rho1 = state1
    else:
        print("State1 is not a vector or density matrix")

    if np.ndim(state2) == 1:
        state2 = toDensity(state2)
        rho2 = np.dot(state2, state2.conj().transpose())
        pure = 1
    elif state2.shape[1] == state2.shape[0]:
        rho2 = state2
    else:
        print("State2 is not a vector or density matrix")

    rho1 = rho1 /np.trace(rho1)
    rho2 = rho2 / np.trace(rho2)

    rho1 = (rho1+rho1.conj().transpose())/2

    if pure:
        val = np.trace(np.dot(rho1, rho2))
    else:
        tmp = sp.linalg.sqrtm(rho1)
        a = np.dot(tmp, np.dot(rho2, tmp))
        val = (np.trace(sp.linalg.sqrtm(a)))**2
    val = np.real(val)

    # when comparing 2 identical pure state, it will get a value larger than 1,
    if val > 1:
        if val - 1 < 0.000001:
            val = 1.0
        else:
            print("Fidelity larger than 1.")

    return val


"""
    concurrence(rhog)
    Desc: Calculates the concurrence of the input state.

    Parameters
    ----------
    rhog : ndarray
        Density Matrix of the desired state.

    Returns
    -------
    val : float
        The calculated concurrence.
    """
def concurrence(rhog):
    if(rhog.shape[0]>2):
        if min(rhog.shape) == 1:
            rhog = np.dot(rhog.conj(), rhog.transpose())

        zz = np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]])
        rr = np.dot(rhog, np.dot(zz, np.dot(rhog.conj(), zz)))
        r = np.linalg.eig(rr)[0]
        # left = np.linalg.inv(right)
        r = np.real(r)

        tmp = np.sort(np.sqrt(r+0j))
        val = np.real(tmp[3]-tmp[2]-tmp[1]-tmp[0])
        val = np.max([val, 0])

        return val
    else:
        return 'NA'


"""
    tangle(rhog)
    Desc: Calculates the tangle of the input state.

    Parameters
    ----------
    rhog : ndarray
        Density Matrix of the desired state. Tangle is calculated by squaring the concurrence.

    Returns
    -------
    val : float
        The calculated tangle.
    """
def tangle(rhog):
    if(rhog.shape[0]>2):
        c = concurrence(rhog)
        val = c ** 2

        return val
    else:
        return 'NA'

"""
    entropy(rhog)
    Desc: Calculates the Von Neumann Entropy of the input state.

    Parameters
    ----------
    rhog : ndarray
        Density Matrix of the desired state.

    Returns
    -------
    val : float
        The calculated entropy.
    """
def entropy(rhog):
    d = np.linalg.eig(rhog)[0]
    e = np.real(d)
    val = 0
    for a in range(len(e)):
        if e[a] > 0:
            val = val-e[a]*np.log2(e[a])

    return val

"""
    linear_entropy(rhog)
    Desc: Calculates the linear entropy of the input state.

    Parameters
    ----------
    rhog : ndarray
        Density Matrix of the desired state.

    Returns
    -------
    val : float
        The calculated linear entropy ranging from 0 to 1/(2^numQubits).
        A value of zero corresponds to a completly pure state.
    """
def linear_entropy(rhog):

    return 1- purity(rhog)


"""
    negativity(rhog)
    Desc: Calculates the negativity of the input state.

    Parameters
    ----------
    rhog : ndarray
        Density Matrix of the desired state.

    Returns
    -------
    val : float
        The calculated negativity.
    """
def negativity(rhog):
    if(rhog.shape[0]>2):
        if min(rhog.shape) == 1:
            rhog = np.dot(rhog, rhog.conj().transpose())

        rho1 = partial_transpose(rhog)
        val = -2*np.min(np.min(np.real(np.linalg.eig(rho1)[0])), 0)

        return val
    else:
        return 'NA'

"""
    purity(rhog)
    Desc: Calculates the purity of the input state.

    Parameters
    ----------
    rhog : ndarray
        Density Matrix of the desired state.

    Returns
    -------
    val : float
        The calculated purity ranging from 1/(2^numQubits) to 1.
        A value of one corresponds to a completly pure state.
    """
def purity(rhog):
    return np.real(np.trace(np.dot(rhog, rhog)))

"""
    partial_transpose(rhog)
    Desc: Returns the partial transpose of the input density matrix.

    DISCLAIMER : In Progress, not checked.

    Parameters
    ----------
    rhog : ndarray
        Input arrays find the partial transpose of.

    Returns
    -------
    rv : ndarray
        Partial transpose of rhog.
    """
def partial_transpose(rhog, n = 0, d = np.nan):
    if min(rhog.shape) == 1:
            rhog = np.dot(rhog, rhog.conj().transpose())

    if d is np.nan:
        n_qubit = int(np.log2(rhog.shape[0]))
        if not n_qubit % 1:
            d = 2*np.ones(n_qubit)
        else:
            print('dimension of rho is incorrect.')

    na = 0
    nb = 0
    nc = 0
    if n < 0:
        na = 1.0
        nb = 1.0
        nc = np.prod(d)
    elif n == 0:
        na = 1.0
        nb = d[n]
        nc = np.prod(d[np.arange(n+1, len(d))])
    elif (n > 0) & (n < len(d)):
        na = np.prod(d[range(n-1)])
        nb = d[n]
        nc = np.prod(d[np.arange(n+1, len(d))])
    elif n == len(d):
        na = np.prod(d[0:n-1])
        nb = d[n]
        nc = 1.0
    elif n > len(d):
        na = np.prod(d)
        nb = 1.0
        nc = 1.0

    if na == 1:
        rv = partial_transpose_helper(rhog, nb)
    # I did't check from here
    else:
        sub_sizes = nb*nc
        y = np.zeros([sub_sizes, sub_sizes, na, na])+0j
        for j in range(sub_sizes):
            for k in range(sub_sizes):
                y[j, k] = rhog[j*sub_sizes:j*sub_sizes+na, k*sub_sizes:k*sub_sizes+na]

        rv = np.zeros([len(rhog), len(rhog)])+0j

        for j in range(na):
            for k in range(na):
                rv[j*nb:j*nb+na, k*nb:k*nb+na] = partial_transpose_helper(y[j, k], nb)

    return rv

"""
    performOperation(psi, g)
    Desc: Performs the operations on the input State.

    Parameters
    ----------
    psi : ndarray
        The input state to do the operation on. Can be a pure state or a density matrix.
    g : ndarray with shape = (num operations, 2^numQubits, 2^numQubits)
        The operations you would like to be done. Can be one operation or an array of operations.

    Returns
    -------
    p : ndarray
        The output state after the operations. Will retain the input form.
    """
def performOperation(psi, g):
    p = psi
    if(len(g.shape) == 3):
        # Multiple operations
        if (len(psi.shape) == 1):
            # ket form
            for i in range(0, len(g)):
                p = ketOperation(p, g[i])
        else:
            # density matrix form
            for i in range(0, len(g)):
                p = densityOperation(p, g[i])
    else:
        # 1 operation
        if (len(psi.shape) == 1):
            # ket form
            p = ketOperation(p, g)
        else:
            # density matrix form
            p = densityOperation(p, g)
    return p


"""
    random_pure_state(N)
    Desc: Returns a random quantum state from a uniform distribution across the space.

    Parameters
    ----------
    N : int
        The dimension of the quantum state

    Returns
    -------
    pure_state : ndarray with shape = (2^numQubits, 2^numQubits)
        The random quantum state.
    """
def random_pure_state(N = 1):
    pure_state = np.zeros(2**N,dtype=complex)
    for x in range(len(pure_state)):
        pure_state[x] = np.random.normal(0, 1)+np.random.normal(0, 1)*1j
    length = np.inner(pure_state.conj(),pure_state)
    pure_state = pure_state/np.sqrt(length)
    length = np.inner(pure_state.conj(), pure_state)
    return pure_state