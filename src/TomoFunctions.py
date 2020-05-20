import scipy as sp
# from numpy.core.defchararray import add
import numpy as np


########################
# TOMOGRAPHY CALCULATE #
########################

#
# def i2array(i, ii, n):
#     nn = np.int(np.ceil((np.log(ii)/np.log(n))))
#     rv = np.zeros(nn)
#     for j in range(nn):
#         rv[j] = i/(n**(nn-j-1))
#         i %= n**(nn-j-1)
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



"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
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


# Converts a density matrix into a lower t matrix
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
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
        #switched order of temp and temp.conj and transpose()
        recurse = np.hsplit(rhog[0:(d-1)], [d-1, d])[0] - np.outer(temp.conj().transpose(), temp)/last_element
    else:
        tm[d-1][0:(d-1)] = np.zeros(d)
        recurse = np.hsplit(rhog[0:(d-1)], [d-1, d])[0]
    for i in range(d-1):
        tm[i][0:(d-1)] = density2tm(recurse)[i][0:(d-1)]

    return tm

# Converts a density matrix into a list of t values
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
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

# Converts a pure state into a density matrix
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def toDensity(psiMat):
    if isinstance(psiMat.size,int):
        return np.outer(psiMat.conj(), psiMat)
    else:
        temp = psiMat[0]
        for j in range(1, len(psiMat)):
            temp = np.kron(temp, psiMat[j])
        return np.outer(temp.conj(), temp)

#
# def one_in(idx, length):
#     val = np.zeros(length)
#     val[idx] = 1
#
#     return val


# Converts a list of t values to an lower t matrix
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
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

# Converts a list of t values to a density matrix
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def t_to_density(t):
    tm = t_matrix(t)
    tm = tm.conj()
    rhog = np.dot(tm.conj().transpose(),tm)

    return rhog


##################
# ERROR ESTIMATE #
##################


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def fidelity(state1, state2):
    pure = 0
    rho1 = state1
    rho2 = state2
    if np.ndim(state1) == 1:
        rho1 = np.dot(state1, state1.conj().transpose())
        pure = 1
    elif state1.shape[1] == state1.shape[0]:
        rho1 = state1
    else:
        print("State1 is not a vector or density matrix")

    if np.ndim(state2) == 1:
        rho2 = np.dot(state2, state2.conj().transpose())
        pure = 1
    elif state2.shape[1] == state2.shape[0]:
        rho2 = state2
    else:
        print("State2 is not a vector or density matrix")

    rho1 /= np.trace(rho1)
    rho2 /= np.trace(rho2)

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
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def concurrence(rhog):
    if(rhog.shape[0]>2):
        if min(rhog.shape) == 1:
            rhog = np.dot(rhog.conj(), rhog.transpose())

        zz = np.array([[0, 0, 0, -1], [0, 0, 1, 0], [0, 1, 0, 0], [-1, 0, 0, 0]])
        rr = np.dot(rhog, np.dot(zz, np.dot(rhog.conj(), zz)))
        r = np.linalg.eig(rr)[0]
        #left = np.linalg.inv(right)
        r = np.real(r)

        tmp = np.sort(np.sqrt(r+0j))
        c = np.real(tmp[3]-tmp[2]-tmp[1]-tmp[0])
        c = np.max([c, 0])

        return c
    else:
        return 0


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def tangle(rhog):
    if(rhog.shape[0]>2>1):
        c = concurrence(rhog)
        t = c ** 2

        return t
    else:
        return 0


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def entanglement(rhog):
    if (rhog.shape[0]>2 > 1):
        t = tangle(rhog)
        x = (1 + np.sqrt(1 - t)) / 2
        if x == 0:
            e = 0
        elif x == 1:
            e = 1
        else:
            e = -x * np.log2(x) - (1 - x) * np.log2(1 - x)

        return e
    else:
        return 0


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def entropy(rhog):
    d = np.linalg.eig(rhog)[0]
    e = np.real(d)
    s = 0
    for a in range(len(e)):
        if e[a] > 0:
            s = s-e[a]*np.log2(e[a])

    return s


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def linear_entropy(rhog):
    if min(rhog.shape) == 1:
        lin_e = 0
    else:
        d = len(rhog)
        lin_e = d * np.real(1-np.trace(np.dot(rhog, rhog)))/(d-1)

    return lin_e


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
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


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def partial_transpose(rhog, n, d=np.nan):
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
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def negativity(rhog):
    if min(rhog.shape) == 1:
        rhog = np.dot(rhog, rhog.conj().transpose())

    rho1 = partial_transpose(rhog, 0)
    val = -2*np.min(np.min(np.real(np.linalg.eig(rho1)[0])), 0)

    return val

#
# def trace_dist(rho1, rho2):
#     #didn't checked, and would not be called in this version.
#     s1 = rho2stokes(rho1)
#     s2 = rho2stokes(rho2)
#     s = s1 - s2
#     val = np.sqrt(np.dot(s.conj().transpose(), s))/2
#
#     return val

################
## Bell State ##
################

####################
## OTHER MEASURES ##
####################


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def purity(rhog):
    return np.real(np.trace(np.dot(rhog, rhog)))


"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def rsquare(rhog):
    return (1+purity(rhog))/2

# performs the operation on the density matrix
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def densityOperation(psi, gate):
    return np.matmul(np.matmul(gate, psi), np.conjugate(np.transpose(gate)))

# performs the operation on the ket state
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def ketOperation(psi, gate):
    return np.matmul(gate, psi)

# Performs the operations on a ket state or density matrix
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def performOperation(psi, g):
    p = psi
    if (len(psi.shape) == 1):
        # ket form
        for i in range(0, len(g)):
            p = ketOperation(p, g[i])
    else:
        # density matrix form
        for i in range(0, len(g)):
            p = densityOperation(p, g[i])
    return p

# calculates the projection values of a onto v
"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
def projVal(v, a):
    # projects a onto v
    return np.dot(a, v) / np.dot(v, v)

