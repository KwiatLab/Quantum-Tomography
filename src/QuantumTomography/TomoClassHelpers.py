from __future__ import print_function
import numpy as np
from .TomoFunctions import *

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


# These are various helper functions used in other main functions



###############
# LIKELYHOODS #
###############

"""
maxlike_fitness(t, coincidences, accidentals, m, prediction)
Desc: Calculates the diffrence between the current predicted state data and the actual data.

Parameters
----------
t : ndarray
    T values of the current predicted state.
coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
    The counts of the tomography.
accidentals : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
    The singles values of the tomography. Used for accidental correction.
m : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
    The measurements of the tomography in density matrix form.
prediction : ndarray
    Predicted counts from the predicted state.
overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
    The relative weights of each measurment. Used for drift correction.
    
Returns
-------
val : float
    value of the optimization function.
"""
def maxlike_fitness(t, coincidences, accidentals, measurements, overall_norms):
    rhog = t_to_density(t,normalize=False)
    prediction = np.zeros_like(coincidences)
    for j in range(len(prediction)):
        prediction[j] = overall_norms[j] * np.real(np.trace(np.dot(measurements[j, :, :], rhog))) + accidentals[j]
        prediction[j] = np.max([prediction[j], 0.01])
    log_like = (prediction - coincidences) / np.sqrt(prediction)
    return np.real(log_like)
"""
maxlike_fitness_hedged(t, coincidences, accidentals, measurements, prediction, bet)
Desc: Calculates the diffrence between the current predicted state data and the actual data using hedged maximum likelihood.

Parameters
----------
t : ndarray
    T values of the current predicted state.
coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
    The counts of the tomography.
accidentals : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
    The singles values of the tomography. Used for accidental correction.
measurements : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
    The measurements of the tomography in density matrix form .
prediction : ndarray
    Predicted counts from the predicted state.
bet : float 0 to .5
    The beta value used.
overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
    The relative weights of each measurment. Used for drift correction.

Returns
-------
val : float
    value of the optimization function.
"""
def maxlike_fitness_hedged(t, coincidences, accidentals, measurements, bet,overall_norms):
    prediction = np.zeros_like(coincidences)
    rhog = t_to_density(t,normalize=False)
    for j in range(len(prediction)):
        prediction[j] = overall_norms[j] * np.real(np.trace(np.dot(measurements[j,:,:], rhog))) + accidentals[j]
        prediction[j] = np.max([prediction[j], 0.01])
    hedge = np.repeat(np.real((bet * np.log(np.linalg.det(rhog))) / len(prediction)), len(prediction))
    val = np.sqrt(np.real((((prediction - coincidences) ** 2) / (2 * prediction)) - hedge) + 1000)
    return np.real(val)

# This are some old functions that have since been reworked.
#
# """
# maxlike_fitness_old(t, coincidences, accidentals, m, prediction)
# Desc: Calculates the diffrence between the current predicted state data and the actual data
#
# Parameters
# ----------
# t : ndarray
#     T values of the current predicted state.
# coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
#     The counts of the tomography.
# accidentals : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
#     The singles values of the tomography. Used for accidental correction.
# m : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
#     The measurements of the tomography in density matrix form.
# prediction : ndarray
#     Predicted counts from the predicted state.
# overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
#     The relative weights of each measurment. Used for drift correction.
#
# Returns
# -------
# val : float
#     value of the optimization function.
# """
# def maxlike_fitness_old(t, coincidences, accidentals, m, prediction,overall_norms=-1):
#     rhog = t_to_density(t,normalize=False)
#     for j in range(len(prediction)):
#         prediction[j] = np.float64(np.real(overall_norms[j] * np.real(np.trace(np.dot(m[:, :, j], rhog))) + accidentals[j]))
#         prediction[j] = np.max([prediction[j], 0.01])
#     val = (prediction - coincidences) / np.sqrt(prediction)
#     val = np.float64(np.real(val))
#     return val
#
#
# """
# maxlike_fitness_hedged_old(t, coincidences, accidentals, m, prediction, bet)
# Desc: Calculates the diffrence between the current predicted state data and the actual data using hedged maximum likelihood.
#
# Parameters
# ----------
# t : ndarray
#     T values of the current predicted state.
# coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
#     The counts of the tomography.
# accidentals : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
#     The singles values of the tomography. Used for accidental correction.
# m : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
#     The measurements of the tomography in density matrix form .
# prediction : ndarray
#     Predicted counts from the predicted state.
# bet : float 0 to .5
#     The beta value used.
# overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
#     The relative weights of each measurment. Used for drift correction.
#
# Returns
# -------
# val : float
#     value of the optimization function.
# """
# def maxlike_fitness_hedged_old(t, coincidences, accidentals, m, prediction, bet,overall_norms=-1):
#     # If overall_norms not given then assume uniform
#     if not isinstance(overall_norms,np.ndarray):
#         overall_norms = np.ones(coincidences.shape[0])
#     elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
#         raise ValueError("Invalid intensities array")
#
#     rhog = t_to_density(t,normalize=False)
#     for j in range(len(prediction)):
#         prediction[j] = overall_norms[j] * np.real(np.trace(np.dot(m[:, :, j], rhog))) + accidentals[j]
#         prediction[j] = np.max([prediction[j], 0.01])
#
#     hedge = np.repeat(np.real((bet * np.log(np.linalg.det(np.mat(rhog)))) / len(prediction)), len(prediction))
#     val = np.sqrt(np.real((((prediction - coincidences) ** 2) / (2 * prediction)) - hedge) + 1000)
#
#     val = np.float64(np.real(val))
#
#     return val

# helper function used in linear tomography
def independent_set(measurements):
    m = measurements[0, :].conj().transpose()
    # may have to switched order of m and measurements, may be wrong but has little effect
    matrix = rho2stokes(np.outer(m, measurements[0, :]))
    max_rank = matrix.shape[0]

    if (measurements.shape[0]) == max_rank:
        s = np.ones([measurements.shape[0], 1])
        return s

    s = np.zeros([measurements.shape[0], 1])
    s[0] = 1
    cur_rank = 1
    for j in np.arange(1, measurements.shape[0], 1):
    # may have to switched order of m and measurements, may be wrong but has little effect
        m = measurements[j, :].conj().transpose()
        sv = rho2stokes(np.outer(m, measurements[j, :]))
        if (np.linalg.matrix_rank(np.concatenate((matrix, sv), axis = 1), tol = 0.001)) > cur_rank:
            matrix = np.concatenate((matrix, sv), axis = 1)
            cur_rank += 1
            s[j] = 1
        else:
            s[j] = 0
        if cur_rank == max_rank:
            break
    return s

# This function converts a list of loglikelihoods to normalized likelihoods. Also returns the index of
# the min loglike
def normalizeLikelihoods(likelihoods):
    nIndex = np.argmin(likelihoods)
    nFactor = likelihoods[nIndex]
    scaled = likelihoods - nFactor
    likelihoods = np.exp(-1 * scaled)
    likelihoods = likelihoods/sum(likelihoods)
    return likelihoods, nIndex,scaled


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

# # Helper function for looping through a multi indexed array
# def multiloop_index(j, lengths):
#     ind = np.zeros(len(lengths))
#     for k in range(len(lengths)-1):
#         sz = np.prod(lengths[np.arange(k+1, len(lengths))])
#         ind[k] = np.fix(j/sz)+1
#         j % = sz
#     ind[len(ind)-1] = j+1
#
#     return ind

def make_positive(rhog_in):
    d, v = np.linalg.eig(rhog_in)
    rhog = np.zeros(rhog_in.shape)
    for j in range(len(d)):
        rhog = rhog + np.abs(d[j])*np.outer(v[:, j], v[:, j].conj().transpose())
    rhog = (rhog + rhog.conj().transpose())/2.0
    return rhog

# helper function for re formatting matrices
def sigma_n(j, nn):
    if j < 0 or j > nn**2-1:
        print('sigma_N: j out of range for SU(N)')

    m = int(np.fix(j/nn))
    n = int(j % nn)
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
    else:  # n = m = N
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

# Evaluates the given function with the given arguments
def fevel(funcname, *args):
    return eval(funcname)(*args)

# Helper function that returns the given properties of the given rho.
# This will not given you fval and intensity. So it is recommended to use getProperties in the TomoClass
def getProperties_helper(errf, rho0):
    return np.array([[errf, fevel(errf, rho0)]for errf in errf], dtype = 'O')


# Helper function that returns the given properties with bounds of the given rhos.
# This will not given you fvals and intensities. So it is recommended to use getProperties in the TomoClass
def getProperties_helper_bounds(errf, rhos):
    n = rhos.shape[0]
    n_fun = len(errf)
    data = np.zeros([n, n_fun], dtype = "O")
    for j in range(n):
        data[j] = getProperties_helper(errf, rhos[j])[:, 1]

    stds = np.zeros(n_fun, dtype = "O")
    means = np.zeros(n_fun, dtype = "O")

    for k in range(n_fun):
        if all(data[:, k] == 'NA'):
            stds[k] = 'NA'
            means[k] = 'NA'
        else:
            if n < 2:
                stds[k] = "NA"
                means[k] = np.mean(data[:, k])
            else:
                stds[k] = np.std(data[:, k], ddof = 1)
                means[k] = np.mean(data[:, k])

    return np.array([np.array(errf), means, stds], dtype = "O").transpose()

# used in getBellSettings
def bellsettings_range_init(rhog, partsize):
    sval = 0
    aval = 0
    apval = 0
    bval = 0
    bpval = 0

    for a in np.linspace(0, np.pi / 2, partsize):
        for ap in np.linspace(a, np.pi / 2, partsize):
            for b in np.linspace(0, np.pi / 2, partsize):
                for bp in np.linspace(b, np.pi / 2, partsize):
                    npp = np.real(np.trace(np.dot(coinmat(a, b), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, b + np.pi / 2), rhog)))
                    e_ab = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, b), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, b + np.pi / 2), rhog)))
                    e_apb = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(a, bp), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, bp + np.pi / 2), rhog)))
                    e_abp = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, bp), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, bp + np.pi / 2), rhog)))
                    e_apbp = 2 * (npp + nmm) - 1

                    s = e_ab + e_abp + e_apb + e_apbp - 2 * np.min([e_ab, e_abp, e_apb, e_apbp])

                    if s > sval:
                        sval = s
                        aval = a
                        apval = ap
                        bval = b
                        bpval = bp

    arange_s = [np.max([aval - ((np.pi / 2) / partsize), 0]), np.min([aval + ((np.pi / 2) / partsize), np.pi / 2])]
    aprange_s = [np.max([apval - ((np.pi / 2) / partsize), 0]),
                 np.min([apval + ((np.pi / 2) / partsize), np.pi / 2])]
    brange_s = [np.max([bval - ((np.pi / 2) / partsize), 0]), np.min([bval + ((np.pi / 2) / partsize), np.pi / 2])]
    bprange_s = [np.max([bpval - ((np.pi / 2) / partsize), 0]),
                 np.min([bpval + ((np.pi / 2) / partsize), np.pi / 2])]

    return [sval, arange_s, brange_s, aprange_s, bprange_s]

# used in getBellSettings
def bellsettings_range(rhog, partsize, arange, brange, aprange, bprange):

    sval = 0
    aval = 0
    apval = 0
    bval = 0
    bpval = 0

    for a in np.linspace(arange[0], arange[1], partsize):
        for ap in np.linspace(aprange[0], aprange[1], partsize):
            for b in np.linspace(brange[0], brange[1], partsize):
                for bp in np.linspace(bprange[0], bprange[1], partsize):
                    npp = np.real(np.trace(np.dot(coinmat(a, b), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, b + np.pi / 2), rhog)))
                    e_ab = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, b), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, b + np.pi / 2), rhog)))
                    e_apb = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(a, bp), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, bp + np.pi / 2), rhog)))
                    e_abp = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, bp), rhog)))
                    nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, bp + np.pi / 2), rhog)))
                    e_apbp = 2 * (npp + nmm) - 1

                    s = e_ab + e_abp + e_apb + e_apbp - 2 * np.min([e_ab, e_abp, e_apb, e_apbp])

                    if s > sval:
                        sval = s
                        aval = a
                        apval = ap
                        bval = b
                        bpval = bp

    arange_s = [np.max([aval - ((arange[1] - arange[0]) / partsize), 0]),
                np.min([aval + ((arange[1] - arange[0]) / partsize), np.pi / 2])]
    aprange_s = [np.max([apval - ((aprange[1] - aprange[0]) / partsize), 0]),
                 np.min([apval + ((aprange[1] - aprange[0]) / partsize), np.pi / 2])]
    brange_s = [np.max([bval - ((brange[1] - brange[0]) / partsize), 0]),
                np.min([bval + ((brange[1] - brange[0]) / partsize), np.pi / 2])]
    bprange_s = [np.max([bpval - ((bprange[1] - bprange[0]) / partsize), 0]),
                 np.min([bpval + ((bprange[1] - bprange[0]) / partsize), np.pi / 2])]

    return [sval, aval, apval, bval, bpval, arange_s, brange_s, aprange_s, bprange_s]

# Helper function that returns the optimal CHSH bell measurement settings given rho.
# It is recommended to use getBellSettings in the TomoClass
def getBellSettings_helper(rhog, partsize_init, partsize, t):

    [s, arange, brange, aprange, bprange] = bellsettings_range_init(rhog, partsize_init)
    a = 0
    b = 0
    ap = 0
    bp = 0

    for j in range(t):
        [s, a, ap, b, bp, arange, brange, aprange, bprange] = \
            bellsettings_range(rhog, partsize, arange, brange, aprange, bprange)

    return np.array([["s", s], ["a", a], ["a'", ap], ["b", b], ["b'", bp]], dtype = "O")


# Helper function that returns the optimal CHSH bell measurement settings with bounds of the given rhos.
# It is recommended to use getBellSettings in the TomoClass
def getBellSettings_helper_bounds(rhop, rho, partsize_init, partsize, t, n):
    belldata = np.zeros([n+1, 5])

    for j in range(n):
        belldata[j] = getBellSettings_helper(rhop[j], partsize_init, partsize, t)[:, 1]
    [bellNames, belldata[-1]] = getBellSettings_helper(rho, partsize_init, partsize, t).transpose()
    bmeans = np.zeros(5)
    berrors = np.zeros(5)

    for m in range(5):
        berrors[m] = np.std(belldata[:, m], ddof = n - 1)
        bmeans[m] = np.mean(belldata[:, m])

    return np.array([bellNames, berrors, bmeans], dtype = "O").transpose()

# Calculates the weighted covariance. Used in bayesian tomography
def weightedcov(samples,weights):
    mean = np.dot(weights,samples)
    diff = samples-mean
    covariance = np.array([np.outer(d, d) for d in diff])
    covariance = np.dot(covariance.T,weights)
    return [mean,covariance]