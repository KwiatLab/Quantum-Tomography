import numpy as np
import scipy
from .TomoClass import Tomography
from .TomoClassHelpers import t_to_density
from .TomoFunctions import toDensity


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""




"""
Returns the default linearly independent states that span the space for 4 measurements of 6 measurements (default).
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
Gets the density matrices of the inputs, outputs, and measurement states from the coincidences, 
input states, and measurement states.
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
gets c matrix from Joe Altepeter's Thesis page 69 eq. 4.52
"""
def get_c_matrix(output_densities, measurement_densities):
    n_measurements = len(measurement_densities.keys())

    # flattening the measurement densities into a single matrix so that the c matrix can be solved for.
    # each column of this matrix is a 2x2 measurement density flattened into its 4 elements
    flattened_measurements = np.zeros((4,n_measurements), dtype=complex)
    for i in range(n_measurements):
        flattened_measurements[:,i] = np.ndarray.flatten(measurement_densities[i])

    c_matrix = np.zeros((n_measurements, n_measurements))

    for i in range(n_measurements):
        current_output_rho_flat = np.ndarray.flatten(output_densities[i])
        c_matrix[i,:] = np.linalg.lstsq(flattened_measurements, current_output_rho_flat)[0]

    return c_matrix

"""
Gets B matrix from Joe Altepeter's thesis page 69 equation 4.53
"""
def get_b_matrix(input_densities, measurement_densities):
    n_measurements = len(measurement_densities.keys())

    # flattening the measurement densities into a single matrix so that the c matrix can be solved for.
    # each column of this matrix is a 2x2 measurement density flattened into its 4 elements
    flattened_measurements = np.zeros((4, n_measurements), dtype=complex)
    for i in range(n_measurements):
        flattened_measurements[:, i] = np.ndarray.flatten(measurement_densities[i])

    b_matrix = np.zeros((4, 4, n_measurements, n_measurements), dtype=complex)
    paulis = get_paulis()

    # equation 4.53 in Joe Altepeter's thesis, solving for the elements of B
    for m in range(4):
        for n in range(4):
            for j in range(n_measurements):
                current_input = paulis[m] @ input_densities[j] @ paulis[n].conj().transpose()
                current_input_flat = np.ndarray.flatten(current_input)

                b_matrix[m, n, j, :] = np.linalg.lstsq(flattened_measurements, current_input_flat)[0]

    return b_matrix

"""
Constructs the chi matrix from the B and C matrices labeled above.
"""
def construct_chi(b_matrix, c_matrix):
    n_measurements = b_matrix.shape[2]

    inv_b = np.linalg.pinv(b_matrix)
    chi_matrix = np.zeros((4,4), dtype=complex)

    for m in range(4):
        for n in range(4):
            for j in range(n_measurements):
                for k in range(n_measurements):
                    chi_matrix[m,n] += inv_b[m,n,j,k] * c_matrix[j,k]

    chi_matrix = chi_matrix / np.trace(chi_matrix)
    return chi_matrix





"""
Returns the pauli matrices for use in Process Tomography
"""
def get_paulis():
    default_paulis = {
        0: np.array([[1, 0], [0, 1]], dtype=complex),
        1: np.array([[0, 1], [1, 0]], dtype=complex),
        2: np.array([[0, -1j], [1j, 0]], dtype=complex),
        3: np.array([[1, 0], [0, -1]], dtype=complex)
    }
    return default_paulis