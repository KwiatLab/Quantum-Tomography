import numpy as np
from . import tomo_math, state_utilities

###############
# LIKELYHOODS #
###############

class Tomography():
    def __init__(self):
        pass

# Helper function for calculating the bell settings
def coinmat(a, b):
    k = np.array([np.cos(a)*np.cos(b), np.cos(a)*np.sin(b), np.sin(a)*np.cos(b), np.sin(a)*np.sin(b)])
    cmat = np.outer(k, k)

    return cmat

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
    rhog = state_utilities.t_to_density(t,normalize=False)
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
    rhog = state_utilities.t_to_density(t,normalize=False)
    for j in range(len(prediction)):
        prediction[j] = overall_norms[j] * np.real(np.trace(np.dot(measurements[j,:,:], rhog))) + accidentals[j]
        prediction[j] = np.max([prediction[j], 0.01])
    hedge = np.repeat(np.real((bet * np.log(np.linalg.det(rhog))) / len(prediction)), len(prediction))
    val = np.sqrt(np.real((((prediction - coincidences) ** 2) / (2 * prediction)) - hedge) + 1000)
    return np.real(val)


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