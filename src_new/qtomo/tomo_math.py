import numpy as np
from . import state_utilities
from random import rand

"""
log_likelyhood(t, coincidences, accidentals, m, prediction)
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
def log_likelyhood(intensity, givenState, coincidences, measurements, accidentals, overall_norms=-1):
    # If overall_norms not given then assume uniform
    if not isinstance(overall_norms, np.ndarray):
        overall_norms = np.ones(coincidences.shape[0])
    elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
        raise ValueError("Invalid intensities array")
    # Convert to densities if given tvals
    if (len(givenState.shape) == 1):
        givenState = state_utilities.t_to_density(givenState)
    # Calculate expected averages
    Averages = np.zeros_like(coincidences, dtype=np.float)
    givenState = intensity*givenState
    for j in range(len(Averages)):
        Averages[j] = overall_norms[j] * np.real(np.trace(np.matmul(measurements[j], givenState))) + \
                      accidentals[j]
        # Avoid dividing by zero for pure states
        if (Averages[j] == 0):
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
    likelihoods = likelihoods/sum(likelihoods)
    return likelihoods, nIndex,scaled

# Calculates the weighted covariance. Used in bayesian tomography
def weightedcov(samples,weights):
    mean = np.dot(weights,samples)
    diff = samples-mean
    covariance = np.array([np.outer(d, d) for d in diff])
    covariance = np.dot(covariance.T,weights)
    return [mean,covariance]


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
    complex = Magn*(np.cos(Phase) +  1j*np.sin(Phase))
    return complex

def complexToPhaser(complex):
    Magn = np.absolute(complex)
    Phase = np.angle(complex)
    if(Phase<0):
        Phase += 2*np.pi
    return np.array([Magn,Phase])