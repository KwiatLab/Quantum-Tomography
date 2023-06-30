from __future__ import print_function
import scipy as sp
import numpy as np
from .TomoFunctionsHelpers import *
import numpy.random as rand
from scipy.linalg import cholesky

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


# # # # # # # # # # # # #
# TOMOGRAPHY CALCULATE  #
# # # # # # # # # # # # #



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
        givenState = t_to_density(givenState)
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
    return np.linalg.cholesky(rhog)

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
    if not isStateVector(psiMat):
        raise ValueError("Invalid input state with shape " + str(psiMat.shape))
    psiMat = np.outer(psiMat, psiMat.conj())
    return psiMat / np.trace(psiMat)

pauli_matrices = np.array([[[1,0],[0,1]],
                           [[0,1],[1,0]],
                           [[0,-1j],[1j,0]],
                           [[1,0],[0,-1]],],dtype=complex)
"""
generalized_pauli_basis(num_qubits)
Desc: Returns a set of pauli matrices.

Parameters
----------
num_qubits : int
    Number of qubits to define the dimensionality of the basis.

Returns
-------
basis : ndarray with shape = (4^numQubits, 2^numQubits, 2^numQubits)
    Pauli Basis
    """
def generalized_pauli_basis(num_qubits):
    if num_qubits == 1:
        return pauli_matrices
    else:
        return np.kron(generalized_pauli_basis(num_qubits-1),pauli_matrices)


"""
get_stokes_parameters(state,basis)
Desc: Given a density or pure state, return the stokes parameters.

Parameters
----------
state : Pure State or Density
    Given State
basis : ndarray with shape = (4^numQubits, 2^numQubits, 2^numQubits)
    Pauli Basis

Returns
-------
stokes : ndarray with length = 4^numQubits
    The stokes parameters for state
    """
def get_stokes_parameters(state,basis=None):
    if basis is None:
        num_qubits = np.log2(state.shape[0])
        num_qubits = int(num_qubits)
        basis = generalized_pauli_basis(num_qubits)
    if isStateVector(state):
        state = toDensity(state)
    stokes = np.array([np.trace(np.matmul(state,b)) for b in basis])
    return stokes

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
    d = int(np.sqrt(len(t)))
    idx = 0
    cur_length = d
    tm = np.zeros([d, d])
    for j in range(int(d)):
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
normalize : bool
    Set to true if you want the density matrix to be normalized by its trace. 
    Default is True.

Returns
-------
rhog : ndarray
    Density Matrix.
    """
def t_to_density(t,normalize=True):
    tm = t_matrix(t)
    rhog = np.dot(tm,tm.T.conj())
    if normalize:
        rhog = rhog / np.trace(rhog)
    return rhog

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
    rho1 = state1
    rho2 = state2
    pure = 0
    if isStateVector(state1):
        rho1 = toDensity(state1)
        pure = 1
    elif state1.shape[1] == state1.shape[0]:
        rho1 = state1
    else:
        raise ValueError("State1 is not a vector or density matrix")

    if isStateVector(state2):
        rho2 = toDensity(state2)
        pure = 1
    elif state2.shape[1] == state2.shape[0]:
        rho2 = state2
    else:
        raise ValueError("State2 is not a vector or density matrix")

    rho1 = rho1 /np.trace(rho1)
    rho2 = rho2 / np.trace(rho2)

    if pure:
        val = np.trace(np.dot(rho1, rho2))
    else:
        tmp = sp.linalg.sqrtm(rho1)
        a = np.dot(tmp, np.dot(rho2, tmp))
        val = (np.trace(sp.linalg.sqrtm(a)))**2
    val = np.real(val)

    # when comparing 2 identical pure state, it will get a value larger than 1,
    if val > 1:
        if val - 1 < 0.00001:
            val = 1.0
        else:
            print("Fidelity larger than 1.")

    return val


# # # # # # # #
# PROPERTIES  #
# # # # # # # #


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
   
Other Properties
 -------------- 
entropy;linear_entropy;negativity;purity;tangle

See Also
 ------ 
err_functions;getProperties
"""
def concurrence(rhog):
    if rhog.shape[0] == 4:
        if isStateVector(rhog):
            rhog = toDensity(rhog)
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

Other Properties
 -------------- 
entropy;linear_entropy;negativity;purity;concurrence

See Also
 ------ 
err_functions;getProperties
"""
def tangle(rhog):
    if rhog.shape[0] == 4:
        if isStateVector(rhog):
            rhog = toDensity(rhog)
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

Other Properties
 -------------- 
concurrence;linear_entropy;negativity;purity;tangle

See Also
 ------ 
err_functions;getProperties
"""
def entropy(rhog):
    if isStateVector(rhog):
        rhog = toDensity(rhog)
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

Other Properties
 -------------- 
entropy;concurrence;negativity;purity;tangle

See Also
 ------ 
err_functions;getProperties
"""
def linear_entropy(rhog):
    if isStateVector(rhog):
        rhog = toDensity(rhog)
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

Other Properties
 -------------- 
entropy;linear_entropy;concurrence;purity;tangle

See Also
 ------ 
err_functions;getProperties
"""
def negativity(rhog):
    if isStateVector(rhog):
        rhog = toDensity(rhog)
    if rhog.shape[0] == 4:
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

Other Properties
 -------------- 
entropy;linear_entropy;negativity;concurrence;tangle

See Also
 ------ 
err_functions;getProperties
"""
def purity(rhog):
    return np.real(np.trace(np.dot(rhog, rhog)))

"""
partial_transpose(rhog)
Desc: Returns the partial transpose of the input density matrix. DISCLAIMER : Tests in progress

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
        for i in range(g.shape[0]):
            p = performOperation(p,g[i])
    else:
        if isStateVector(psi):
            p = ketOperation(p, g)
        else:
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
pure_state : ndarray with length = 2^N
    The random quantum state.
"""
def random_pure_state(N = 1):
    pure_state = np.zeros(2**N,dtype=complex)
    for x in range(len(pure_state)):
        pure_state[x] = rand.normal(0, 1)+rand.normal(0, 1)*1j
    length = np.inner(pure_state.conj(),pure_state)
    pure_state = pure_state/np.sqrt(length)
    length = np.inner(pure_state.conj(), pure_state)
    return pure_state


"""
random_density_state(N)
Desc: Returns a random quantum density from an approximate uniform distribution across the space.

Parameters
----------
N :int
    The dimension of the quantum state

Returns
-------
density : ndarray with shape = (2^N, 2^N)
    The random quantum state.
"""
def random_density_state(N=1):
    density = random_ginibre(2 ** N)
    density = np.matmul(density, density.T.conj())
    if np.trace(density).real<10**-6:
        # If trace is close to 0 then resample
        return random_density_state(N)
    else:
        density = density / np.trace(density)
    return density

"""
random_bell_state(N)
Desc: Randomly returns one of the 4 bell state. 
For 1 qubits on of the standard basis states is returned. For states with dimension
greater then 2 the Greenberger-Horne-Zeilingerstate is returned with a random phase.

Parameters
----------
N :int
    The dimension of the quantum state

Returns
-------
pure_state : ndarray with shape = (2^N, 2^N)
    The random quantum state.
"""
def random_bell_state(N=2):
    if N ==1:
        whichState = rand.randint(0,6)
        if whichState == 0:
            return np.array([1,0],dtype=complex)
        elif whichState == 1:
            return np.array([0,1],dtype=complex)
        elif whichState == 2:
            return 1/np.sqrt(2)*np.array([1,1],dtype=complex)
        elif whichState == 3:
            return 1/np.sqrt(2)*np.array([1,-1],dtype=complex)
        elif whichState == 4:
            return 1/np.sqrt(2)*np.array([1,1j],dtype=complex)
        elif whichState == 5:
            return 1/np.sqrt(2)*np.array([1,-1j],dtype=complex)
    if N ==2:
        whichState = rand.randint(0, 4)
        if whichState == 0:
            return 1 / np.sqrt(2) * np.array([1,0,0,1], dtype=complex)
        elif whichState == 1:
            return 1 / np.sqrt(2) * np.array([1,0,0,-1], dtype=complex)
        elif whichState == 2:
            return 1 / np.sqrt(2) * np.array([0,1,1,0], dtype=complex)
        elif whichState == 3:
            return 1 / np.sqrt(2) * np.array([0,1,-1,0], dtype=complex)
    GHZ = np.zeros(2**N,dtype=complex)
    GHZ[0] = 1
    GHZ[-1] = 1
    return 1 / np.sqrt(2)* GHZ


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

"""
densityOperation(D)
Desc: Performs the operation on the density matrix

Parameters
----------
psi : ndarray with shape = (2^nQubits, 2^nQubits)
    The state to apply the operation to in density form.
gate : ndarray with shape = (2^nQubits, 2^nQubits)
    The operation to apply to the state.

Returns
-------
psi : ndarray with shape = (2^nQubits, 2^nQubits)
    The state with the operation applied.
"""
def densityOperation(psi, gate):
    return np.matmul(np.matmul(gate, psi), np.conjugate(np.transpose(gate)))

"""
ketOperation(D)
Desc: Performs the operation on the pure state.

Parameters
----------
psi : 1darray with length = 2^nQubits
    The state to apply the operation to in ket form.
gate : ndarray with shape = (2^nQubits, 2^nQubits)
    The operation to apply to the state.

Returns
-------
psi : ndarray with shape = (2^nQubits, 2^nQubits)
    The state with the operation applied.
"""
def ketOperation(psi, gate):
    return np.matmul(gate, psi)

"""
quarterWavePlate(D)
Desc: returns the quantum gate associated with a quarter wave plate.

Parameters
----------
theta : float
    The angle of the wave plate with respect to horizontal.
"""
def quarterWavePlate(theta):
    return np.array([
        [np.cos(theta)**2+1j*np.sin(theta)**2   ,   (1-1j)*np.cos(theta)*np.sin(theta)  ],
        [(1-1j)*np.cos(theta)*np.sin(theta)     ,   np.sin(theta)**2+1j*np.cos(theta)**2]
    ],dtype=complex)

"""
halfWavePlate(D)
Desc: returns the quantum gate associated with a half wave plate.

Parameters
----------
theta : float
    The angle of the wave plate with respect to horizontal.
"""
def halfWavePlate(theta):
    return np.array([
        [np.cos(theta)**2-np.sin(theta)**2      ,   2*np.cos(theta)*np.sin(theta)  ],
        [2*np.cos(theta)*np.sin(theta)          ,   np.sin(theta)**2-np.cos(theta)**2]
    ],dtype=complex)

"""
getWavePlateBasis(theta_qwp,theta_hwp,flipPBS)
Desc: Given the angles for the QWP and HWP plate find the measurement basis. PBS is assumed to transmit 
Horizontally polarized light and reflect vertical. This function does not take into account crosstalk.

Parameters
----------
theta_qwp : string
    The angle with respect to horizontal for the quarter wave plate.
theta_hwp : string
    The angle with respect to horizontal for the quarter wave plate.    
flipPBS: bool
        Set this to true to assume the PBS transmits V and reflects H

Returns
-------
basis : ndarray with shape (2,2)
    Top row : State the original state was projected onto given that it transmitted through the PBS,
    Bottom row : State the original state was projected onto given that it reflected off the PBS
"""
def getWavePlateBasis(theta_qwp,theta_hwp,flipPBS=False):
    basis = np.eye(2,dtype=complex)
    if(flipPBS):
        basis = np.fliplr(basis)
    basis = np.matmul(halfWavePlate(theta_hwp).T.conj(),basis)
    basis = np.matmul(quarterWavePlate(theta_qwp).T.conj(),basis)
    # return the transpose so that basis[0] returns state 1
    # no need to do conj, its just a reordering of the numbers
    return basis.T

"""
removeGlobalPhase(D)
Desc: Factors out the global phase of the given state by dividing the entire state by the phase of the first component.

Parameters
----------
pure_state : 1darray with length = 2^nQubits
    The state in ket form.

Returns
-------
psi : ndarray with shape = (2^nQubits, 2^nQubits)
    The state with global phase factored out.
"""
def removeGlobalPhase(pure_state):
    # Normalize
    norm = np.dot(pure_state.conj(), pure_state)
    pure_state = pure_state / norm
    # get phase of first number
    [magn,phase] = complexToPhaser(pure_state[0])
    complexPhaser = phaserToComplex(np.array([1,phase]))
    # divide state by first phase
    pure_state = pure_state/complexPhaser
    return pure_state