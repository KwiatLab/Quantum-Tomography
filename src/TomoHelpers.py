#These are various helper functions used in other main functions

import numpy as np
# helper function used in linear tomography
def independent_set(measurements):
    m = measurements[0, :].conj().transpose()
    #may have to switched order of m and measurements, may be wrong but has little effect
    matrix = rho2stokes(np.outer(m,measurements[0, :]))
    max_rank = matrix.shape[0]

    if (measurements.shape[0]) == max_rank:
        s = np.ones([measurements.shape[0], 1])
        return s

    s = np.zeros([measurements.shape[0], 1])
    s[0] = 1
    cur_rank = 1
    for j in np.arange(1, measurements.shape[0], 1):
    #may have to switched order of m and measurements, may be wrong but has little effect
        m = measurements[j, :].conj().transpose()
        sv = rho2stokes(np.outer(m,measurements[j, :]))
        if (np.linalg.matrix_rank(np.concatenate((matrix, sv), axis=1), tol=0.001)) > cur_rank:
            matrix = np.concatenate((matrix, sv), axis=1)
            cur_rank += 1
            s[j] = 1
        else:
            s[j] = 0
        if cur_rank == max_rank:
            break

    return s

# helper function that formats the projects into a single matrix
def b_matrix(projectors):
    dim_m = projectors.shape[1]
    dim_b = dim_m**2
    tmp = np.zeros([dim_b, dim_b])+0j
    for i in range(dim_b):
        for j in range(dim_b):
            tmp[i][j] = np.inner(projectors[i], np.inner(sigma_n(j, dim_m), projectors[i].conj().transpose()))
    b = tmp.transpose()

    return b


# helper function that formats the measurements into a single matrix
def m_matrix(mu, projectors, b_inv):
    dim_m = projectors.shape[1]
    dim_b = dim_m**2

    tmp = np.zeros([dim_m, dim_m])
    for j in range(dim_b):
        tmp = tmp + b_inv[mu][j]*sigma_n(j, dim_m)
    m = tmp

    return m

# Helper function for looping through a multi indexed array
def multiloop_index(j, lengths):
    ind = np.zeros(len(lengths))
    for k in range(len(lengths)-1):
        sz = np.prod(lengths[np.arange(k+1, len(lengths))])
        ind[k] = np.fix(j/sz)+1
        j %= sz
    ind[len(ind)-1] = j+1

    return ind


#helper function for re formatting matrices
def sigma_n(j, nn):
    if j < 0 or j > nn**2-1:
        print('sigma_N: j out of range for SU(N)')

    m = np.int(np.fix(j/nn))
    n = np.int(j % nn)
    tmp1 = np.zeros([nn, 1])
    tmp2 = np.zeros([nn, 1])
    tmp1[m] = 1
    tmp2[n] = 1

    if m < n:
        matrix = (np.outer(tmp1, tmp2.conj().transpose())+np.outer(tmp2, tmp1.conj().transpose()))*np.sqrt(nn/2.0)
    elif m > n:
        matrix = 1j*(np.outer(tmp1, tmp2.conj().transpose())-np.outer(tmp2, tmp1.conj().transpose()))*np.sqrt(nn/2.0)
    elif (m+1) < nn:
        z = np.zeros(nn)
        for i in range(m+1):
            z[i] = 1
        matrix = -(np.sqrt(nn/((m+1.0)**2+m+1.0)))*np.diag(z)
        matrix[m+1, m+1] = (m+1.0)*(np.sqrt(nn/((m+1.0)**2+m+1.0)))
    else:  # n=m=N
        matrix = np.identity(nn)

    return matrix

# Helper function for independent_set function
def rho2stokes(rhog):
    if rhog.ndim == 1:
        rhog = np.outer(rhog, rhog.conj().transpose())

    d = len(rhog)
    n = d**2

    ss = np.zeros([n, 1])+0j
    for j in range(n):
        ss[j] = np.trace(np.inner(rhog, sigma_n(j, d)))

    return ss

# Helper function for calculating the bell settings
def coinmat(a, b):
    k = np.array([np.cos(a)*np.cos(b), np.cos(a)*np.sin(b), np.sin(a)*np.cos(b), np.sin(a)*np.sin(b)])
    cmat = np.outer(k, k)

    return cmat