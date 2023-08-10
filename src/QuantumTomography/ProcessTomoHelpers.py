import numpy as np
import scipy
from .TomoClass import Tomography
from .TomoClassHelpers import t_to_density
from .TomoFunctions import toDensity, generalized_pauli_basis


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


"""
get_default_states(num_measurements=6)

Desc: Returns the default linearly independent states that span the space for 4 measurements of 6 measurements (default).

Parameters
----------
num_measurements: int
    Number of measurement states (either 4 or 6).

Returns
-------
default_states: ndarray with shape (num_measurements, 2)
    The default linearly independent states that span the space.
"""
def get_default_states(num_measurements=6):
    # H, V, D, A, R, L
    default_states_6 = np.array([[1,0],[0,1],[0.7071,0.7071],[0.7071,-0.7071],[0.7071,0.7071j],[0.7071,-0.7071j]], dtype=complex)
    # H, V, D, R
    default_states_4 = np.array([[1,0],[0,1],[0.7071,0.7071],[0.7071,0.7071j]], dtype=complex)

    if num_measurements == 6:
        return default_states_6
    elif num_measurements == 4:
        return default_states_4
    else:
        raise ValueError("Default States are only defined for 6 and 4 measurements")

"""
get_process_densities(coincidences, input_states, measurement_states)

Desc: Gets the density matrices of the inputs, outputs, and measurement states from the coincidences, 
input states, and measurement states.

Parameters
----------
coincidences : ndarray with shape = (number of measurements, number of measurements)
    Each column corresponds to an input state, and each row corresponds to a measurement state.
input_states : ndarray with shape = (number of measurements, 2) or (number of measurements, number of measurements)
    Each row is a pure state used as an input to the process.
measurement_states : ndarray with shape = (number of measurements, 2)
    Each row is a pure state used to measure the state after the process.   

Returns
-------
input_densities : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
    The states that are input into the quantum process in density matrix form.
measurement_densities: ndarray with shape = (number of measurements, 2, 2)
    The measurements of the tomography in density matrix form.
output_densities : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
    The states that are output from the quantum process in density matrix form.
"""
def get_process_densities(coincidences, input_states, measurement_states):
    input_densities = {}
    measurement_densities = {}
    output_densities = {}

    n_measurements = coincidences.shape[0]

    if input_states.shape == (n_measurements, 2):
        for i in range(n_measurements):
            input_densities[i] = toDensity(input_states[i,:])
    elif input_states.shape == (n_measurements, n_measurements):
        for i in range(n_measurements):
            #if the input states were given as counts, the input states are the columns, and the measurement states the rows
            #so you take the state tomography of each column to give you the rho of each input state
            input_densities[i] = Tomography(1).StateTomography(measurement_states, input_states[:,i])
    else:
        raise ValueError("input states must be mx2 matrix of pure states or mxm matrix of measurement counts")

    for i in range(n_measurements):
        measurement_densities[i] = toDensity(measurement_states[i,:])

        output_densities[i] = Tomography(1).StateTomography(measurement_states, coincidences[:,i])[0]

    return input_densities, measurement_densities, output_densities


"""
get_c_matrix(output_densities, measurement_densities)

Desc: 

Calculates the c matrix from the equation 

\varepsilon(\hat{\rho}_j) = \sum_k c_{jk}\hat{\rho}_k 

by least-squares approximation.

See also Joe Altpeter's thesis (http://research.physics.illinois.edu/QI/Photonics/theses/altepeter-thesis.pdf) Page 69, eq. 4.52.

Parameters
----------
// TODO: double check these shapes
output_densities : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
    The states that are output from the quantum process in density matrix form.

measurement_densities: ndarray with shape = (number of measurements, 2, 2)
    The measurements of the tomography in density matrix form.        

Returns
-------
// TODO: Figure out the shape
c_matrix : ndarray with shape = 
    The approximated c matrix
"""
def get_c_matrix(output_densities, measurement_densities, num_qubits = 1):

    n_measurements = measurement_densities.shape[0]

    # flattening the measurement densities into a single matrix so that the c matrix can be solved for.
    # each column of this matrix is a 2x2 measurement density flattened into its 4 elements
    flattened_measurements = np.zeros((measurement_densities.shape[1]**2,n_measurements), dtype=complex)
    for i in range(n_measurements):
        flattened_measurements[:,i] = measurement_densities[i].flatten()

    c_matrix = np.zeros((n_measurements, n_measurements), dtype=complex)

    for i in range(n_measurements):
        current_output_rho_flat = output_densities[i].flatten()
        c_matrix[i,:] = np.linalg.lstsq(flattened_measurements, current_output_rho_flat, rcond=None)[0]

    return c_matrix


"""
get_b_matrix(input_densities, measurement_densities)

Desc: 

Calculates the beta matrix from the equation
\tilde{E}_m \hat{\rho}_j \tilde{E}_n^\dagger = \sum_k \beta_{jk}^{mn}\hat{\rho}_k
by least-squares approximation.

See also Joe Altpeter's thesis (http://research.physics.illinois.edu/QI/Photonics/theses/altepeter-thesis.pdf) Page 69, eq. 4.53.

Parameters
----------
// TODO: double check these shapes
input_densities : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
    The states that are input into the quantum process in density matrix form.

measurement_densities: ndarray with shape = (number of measurements, 2, 2)
    The measurements of the tomography in density matrix form.

Returns
-------
b_matrix : ndarray with shape = (4,4, number of measurements, number of measurements)
    The approximated beta matrix
"""
def get_b_matrix(input_densities, measurement_densities, num_qubits=1):

    n_measurements = measurement_densities.shape[0]

    # flattening the measurement densities into a single matrix so that the c matrix can be solved for.
    # each column of this matrix is a 2x2 measurement density flattened into its 4 elements
    flattened_measurements = np.zeros((measurement_densities.shape[1]**2,n_measurements), dtype=complex)
    for i in range(n_measurements):
        flattened_measurements[:,i] =measurement_densities[i].flatten()

    b_matrix = np.zeros((4**num_qubits, 4**num_qubits, n_measurements, n_measurements), dtype=complex)
    paulis = generalized_pauli_basis(num_qubits) 
    # equation 4.53 in Joe Altepeter's thesis, solving for the elements of B
    for m in range(4**num_qubits):
        for n in range(4**num_qubits):
            for j in range(n_measurements):
                current_input = paulis[m] @ input_densities[j] @ paulis[n].conj().transpose()
                current_input_flat = current_input.flatten()

                b_matrix[m, n, j, :] = np.linalg.lstsq(flattened_measurements, current_input_flat, rcond=None)[0]

    return b_matrix

"""
construct_chi(b_matrix, c_matrix)

Desc: 

Constructs the \chi matrix from the beta and c matrices given by get_b_matrix() and get_c_matrix().

Each element \chi_{mn} is given by

\chi_{mn} = \sum_{jk} (\beta^-1)_{jk}^{mn} c_{jk}

See also Joe Altpeter's thesis (http://research.physics.illinois.edu/QI/Photonics/theses/altepeter-thesis.pdf) Page 70, eq. 4.55.

Parameters
----------
b_matrix : ndarray with shape = (4,4, number of measurements, number of measurements)
    The approximated beta matrix

// TODO: Figure out shape
c_matrix: ndarray with shape = 
    The approximated c matrix
Returns
-------
chi_matrix : ndarray with shape = (4,4, number of measurements, number of measurements)
    The constructed chi matrix
"""
def construct_chi(b_matrix, c_matrix, num_qubits):
    n_measurements = b_matrix.shape[2]

    inv_b = np.linalg.pinv(b_matrix)
    chi_matrix = np.zeros((4**num_qubits,4**num_qubits), dtype=complex)

    for m in range(4**num_qubits):
        for n in range(4**num_qubits):
            for j in range(n_measurements):
                for k in range(n_measurements):
                    chi_matrix[m,n] += inv_b[m,n,j,k] * c_matrix[j,k]

    chi_matrix = chi_matrix / np.trace(chi_matrix)
    return chi_matrix

"""
get_paulis()

Desc: Returns the Pauli matrices for use in Process Tomography

Returns
-------
default_paulis : dictionary with keys [0,1,2,3] and values as ndarray with shape = (2,2)
    A dictionary containing the default Pauli matrices.
"""

def get_paulis():
    default_paulis = {
        0: np.array([[1, 0], [0, 1]], dtype=complex),
        1: np.array([[0, 1], [1, 0]], dtype=complex),
        2: np.array([[0, -1j], [1j, 0]], dtype=complex),
        3: np.array([[1, 0], [0, -1]], dtype=complex)
    }
    return default_paulis