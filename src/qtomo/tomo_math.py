import numpy as np

from . import state_utilities

"""
log_likelihood(t, coincidences, accidentals, m, prediction)
Desc: This is the log likelyhood function. It used in bayesian tomography.

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
overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit. (optional)
    The relative weights of each measurment. Used for drift correction.
    
Returns
-------
val : float
    value of the optimization function.
"""


def log_likelihood(intensity, givenState, coincidences, measurements, accidentals, overall_norms=-1):
    # If overall_norms not given then assume uniform
    if not isinstance(overall_norms, np.ndarray):
        overall_norms = np.ones(coincidences.shape[0])
    elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
        raise ValueError("Invalid intensities array")
    # Convert to densities if given tvals
    if len(givenState.shape) == 1:
        givenState = state_utilities.t_to_density(givenState)
    # Calculate expected averages
    Averages = np.zeros_like(coincidences, dtype=float)
    givenState = intensity * givenState
    for j in range(len(Averages)):
        Averages[j] = overall_norms[j] * np.real(np.trace(np.matmul(measurements[j], givenState))) + accidentals[j]
        # Avoid dividing by zero for pure states
        if Averages[j] == 0:
            Averages[j] == 1
    val = (Averages - coincidences) ** 2 / (2 * Averages)
    return sum(val)


# This function converts a list of loglikelihoods to normalized likelihoods. Also returns the index of
# the min loglike
def normalizeLikelihoods(likelihoods):
    nIndex = np.argmin(likelihoods)
    nFactor = likelihoods[nIndex]
    scaled = likelihoods - nFactor
    likelihoods = np.exp(-1 * scaled)
    likelihoods = likelihoods / sum(likelihoods)
    return likelihoods, nIndex, scaled


# Calculates the weighted covariance. Used in bayesian tomography
def weightedcov(samples, weights):
    mean = np.dot(weights, samples)
    diff = samples - mean
    covariance = np.array([np.outer(d, d) for d in diff])
    covariance = np.dot(covariance.T, weights)
    return [mean, covariance]


def ginibre(rng: np.random.Generator, d: int = 2):
    """Return a dxd matrix with entries from the Ginibre ensemble.

    Matrix entries are of the form a+ib | a,b iid. Norm(0,1)
    """
    a = rng.normal(size=(d, d))
    b = rng.normal(size=(d, d))
    return a + b * 1j
