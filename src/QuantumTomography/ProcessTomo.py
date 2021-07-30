import numpy as np
import scipy
from .ProcessTomoHelpers import *
from .TomoClassHelpers import make_positive
from .TomoFunctions import density2t, t_to_density

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""


"""
ProcessTomography(input_states, measurement_states, coincidences)
Desc: This function calculates the chi matrix most likely to have produced the observed counts

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
chi_matrix : ndarray with shape = (4,4)
    Represents the chi matrix that best characterizes the process.
"""
def ProcessTomography(coincidences, input_states=get_default_states(6), measurement_states=get_default_states(6)):

    # using tomography to calculate the density matrices after the process, and converting the input pure states (or coincidences) and measurement pure states to density matrices
    input_densities, measurement_densities, output_densities = get_process_densities(coincidences,
                                                                                           input_states,
                                                                                           measurement_states)

    # construct the starting chi matrix to be input into the MLE process tomography
    chi_matrix = StandardProcessTomography(input_densities, measurement_densities, output_densities)

    # makes the chi matrix positive definite so that the cholesky decomp can be used in MLE
    chi_matrix = make_positive(chi_matrix)

    # performs the MLE on the chi matrix
    chi_matrix = MLEProcessTomography(chi_matrix, input_densities, measurement_densities, coincidences)

    return chi_matrix


"""
MLEProcessTomography(chi_matrix), coincidences, input_densities, measurement_densities)
Desc: Uses Maximum Likelihood Estimation to find the chi matrix most likely to have produced the observed measurements.

Parameters
----------
chi_matrix : ndarray with shape = (4, 4)
    The chi matrix to start the MLE with. Must be Hermitian.
coincidences : ndarray with shape = (number of measurements, number of inputs)
    A matrix holding the measured coincidences, with the inputs in each column and different measurement states in each row.
input_densities : dictionary with keys 0-number of inputs
    Dictionary that holds the density matrices for each input state. Order must be the same as in the columns of coincidence matrix.
measurement_densities : dictionary with keys 0-number of measurements
    Dictionary that holds the density matrices for each measurement state. Order must be the same as in rows of coincidence matrix.

"""
def MLEProcessTomography(chi_matrix, coincidences, input_densities, measurement_densities):
    #converting the chi matrix to its t values so that there are less parameters for optimization
    chi_matrix_t_vals = density2t(chi_matrix)
    chi_matrix_t_vals += 0.0001

    # get the final t values from scipy.optimize.least sq on the process_count_differences function
    final_chi_t_vals = scipy.optimize.leastsq(process_count_differences, chi_matrix_t_vals, args=(input_densities, measurement_densities, coincidences))[0]

    # converting the t values back into a chi matrix and then normalizing
    chi_matrix = t_to_density(final_chi_t_vals)
    chi_matrix = chi_matrix / np.trace(chi_matrix)

    return chi_matrix


"""
StandardProcessTomography(input_densities, measurement_densities, output_densities)
Desc: Uses Standard Process Tomography (see Joe Altepeter's thesis pg. 69) to calculate a likely chi matrix for a given process.

Parameters
----------
input_densities : dictionary with number of elements = number of input states
    A dictionary (with integer keys starting at 0), with each value representing a density matrix for an input.
measurement_densities : dictionary with number of elements = number of measurements
    A dictionary (again with integer keys starting at 0) with each value representing the density matrix of a measurement state
output_densities : dictionary with number of elements = number of input states
    A dictionary in which each value represents the output for the corresponding input in the input_densities dictionary.
    
Returns
-------
chi_matrix : ndarray with shape = (4, 4)
    The chi matrix that characterizes the process.
"""
def StandardProcessTomography(input_densities, measurement_densities, output_densities):
    # get the number of measurements
    n_measurements = len(measurement_densities.keys())

    # c matrix definted in eq. 4.52 of Joe Altepeter's thesis page 69
    c_matrix = get_c_matrix(output_densities, measurement_densities)

    # B matrix defined in eq. 4.53 of Joe Altepeter's thesis page 69
    b_matrix = get_b_matrix(input_densities, measurement_densities)

    # constructing the chi matrix from the B and c matrices in eq. 4.55 pg. 60 of Joe Altepeter's thesis
    chi_matrix = construct_chi(b_matrix, c_matrix)

    return chi_matrix


"""
post_process_ensity(chi_matrix, rho)
Desc: Calculates the post-process density matrix of a given input rho through a process characterized by the chi matrix

Parameters
----------
chi_matrix : ndarray with shape = (4, 4)
    The chi matrix that characterizes the process to send the input through.
rho : ndarray with shape = (2, 2)
    Density matrix representing the input state to the process.
    
Returns
-------
output_rho : ndarray with shape = (4, 4)
    The density matrix output from the process.
"""
def post_process_density(chi_matrix, rho):
    output_rho = np.zeros_like(rho)
    paulis = get_paulis()

    for m in range(len(paulis.keys())):
        for n in range(len(paulis.keys())):
            output_rho += chi_matrix[m, n] * (paulis[m] @ rho @ paulis[n].conj().transpose())

    return output_rho