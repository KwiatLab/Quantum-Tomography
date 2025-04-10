import numpy as np
import typing
import scipy as sp
import sympy as sym
import json
from sympy import init_printing
import argparse

init_printing()
np.set_printoptions(formatter={"float": lambda x: f"{x:10.4g}"}, suppress=True)

"""
Created on: 2025-04-09
    Authors: Timur Javid
    Kwiat Quantum Information Group

Quantum Process Tomography

This script performs quantum state and process tomography using optimization techniques. 
It reconstructs density matrices and process matrices from experimental measurement data.

User Defined Constants:
- PAULI_I, PAULI_X, PAULI_Y, PAULI_Z: Pauli matrices. TODO: Make this able to take a generic basis
- PAULIS: List of Pauli matrices.
- NUM_QUBITS: Number of qubits (set to 1). TODO: User friendly way to change this.
- BASIN_HOPPING_N_ITER_STATE: Number of iterations for basin hopping in state tomography.
- BASIN_HOPPING_N_ITER_PROCESS: Number of iterations for basin hopping in process tomography.

Main Execution:
- Parses command-line arguments to load input data.
- Prepares input and output counts and measurement densities.
- Performs process tomography to reconstruct the chi matrix.
- Saves the results to a file.
"""

PAULI_I = np.eye(2, dtype=complex)
PAULI_X = np.array([[0, 1], [1, 0]], dtype=complex)
PAULI_Y = np.array([[0, -1.0j], [1.0j, 0]], dtype=complex)
PAULI_Z = np.array([[1, 0], [0, -1]], dtype=complex)

PAULIS = [PAULI_I, PAULI_X, PAULI_Y, PAULI_Z]

NUM_QUBITS = 1

BASIN_HOPPING_N_ITER_STATE = 10
BASIN_HOPPING_N_ITER_PROCESS = 20


def test_density(params: typing.Tuple[float, float, float, float]) -> np.ndarray:
    """
    Constructs a density matrix based on the given real parameters.

    Parameters:
        params (tuple): (i1, i2, i3, i4) 4 real values that
                        parameterize a Hermitian matrix.
    Returns:
        numpy.ndarray: Density matrix constructed from lower
                       triangular parameterization.
    """
    i1, i2, i3, i4 = params
    L = np.array([[i1, 0], [i3 + 1.0j * i4, i2]])
    return L @ L.conj().T


# Output the expected counts for a given test density and input state
def expected_counts(state, density):
    return np.trace(state @ density)


def expected_counts_vec(measurement_densities: np.ndarray, rho: np.ndarray) -> np.ndarray:
    """
    Gives you the expected counts when performing the measurements given
    on a density matrix rho

    Parameters:
        measurement_densities (numpy.ndarray): Measurement density matrices for the measurements you performed.
        rho (numpy.ndarray): Density matrix to be measured.
    Returns:
        numpy.ndarray: Density matrix constructed from lower
                       triangular parameterization.
    """
    total_measurement = measurement_densities.reshape(
        measurement_densities.shape[0], (2**NUM_QUBITS) ** 2
    )
    return np.real(
        total_measurement @ rho.flatten(order="F")
    )  # Fortran order is column order flattening, need it to preserve the multiplication we're doing


def state_cost(
    params: np.ndarray, state_index: int, counts: np.ndarray, measurement_densities: np.ndarray
) -> float:
    """
    Calculates the cost function for state tomography

    Parameters:
        params (numpy.ndarray): Real parameters for the density matrix.
        state_index (int): The index of the state we're measuring w.r.t. the counts matrix.
        measurement_densities (numpy.ndarray): Measurement density matrices for the measurements you performed..
    Returns:
        float: Cost function value, log likelihood of a Poissonian distribution.
               Calculates (expected - actual)^2/(2*expected) for each measurement and sums
    """
    density = test_density(params)
    expected = expected_counts_vec(measurement_densities, density)

    diff_vec = expected - counts[state_index].flatten(
        order="F"
    )  # Need to match the flatten ordering from expected_counts_vec()
    expected_count_normalization = np.diag(1 / (2 * expected))
    return np.abs(diff_vec.T @ expected_count_normalization @ diff_vec)


def normalize_and_chop(tomo_results: np.ndarray) -> np.ndarray:
    """
    Normalizes tomography results and removes small numerical artifacts.

    Parameters:
        tomo_results (numpy.ndarray): Tomography results (assuming it's a list of density matrices).
    Returns:
        numpy.ndarray: Normalized and chopped density matrices.
    """
    traces = np.trace(tomo_results, axis1=1, axis2=2)
    res = np.divide(tomo_results, traces[:, np.newaxis, np.newaxis])
    real_part = np.real(res)
    imag_part = np.imag(res)
    chopped_real = np.where(np.abs(real_part) < 1e-4, 0, real_part)
    chopped_imag = np.where(np.abs(imag_part) < 1e-4, 0, imag_part)
    return chopped_real + 1.0j * chopped_imag


def calculate_fidelity(input: np.ndarray, reference: np.ndarray) -> float:
    """
    Calculates fidelity between two density matrices.

    Parameters:
        input (numpy.ndarray): Input density matrix.
        reference (numpy.ndarray): Reference density matrix.
    Returns:
        float: Fidelity w.r.t. the reference density matrix.
    """
    sqrt_actual = sp.linalg.fractional_matrix_power(input, 0.5)
    multiplied = sqrt_actual @ reference @ sqrt_actual
    sqrt_mult = sp.linalg.fractional_matrix_power(multiplied, 0.5)
    fidelity = np.real(np.trace(sqrt_mult) ** 2)
    return fidelity


# TODO: Make this generic for any number of qubits
def normalize_counts(counts: np.ndarray) -> np.ndarray:
    """
    Normalizes the counts w.r.t. the measurements. Ex: HH/(HH+HV), HV/(HH+HV), etc.

    Parameters:
        counts (numpy.ndarray): Counts matrix.
    Returns:
        numpy.ndarray: Normalized counts matrix.
    """
    counts = np.array(counts)
    for i in range(0, counts.shape[0], 2):
        normalization_factor = counts[:, i] + counts[:, i + 1]
        counts[:, i] = counts[:, i] / normalization_factor
        counts[:, i + 1] = counts[:, i + 1] / normalization_factor
    return counts


def state_tomography(counts, measurement_densities):
    """
    State tomography from measured counts.

    Parameters:
        counts (numpy.ndarray): Measured counts.
        measurement_densities (numpy.ndarray): Measurement density matrices for the measurements you performed.
    Returns:
        np.ndarray: Tomography results (density matrices).
        np.ndarray: Array of purities (float values) of the reconstructed density matrices.
        float: Final cost function value.
    """
    normalized_counts = normalize_counts(counts)

    tomography_results = []
    for i in range(len(measurement_densities)):
        ## Non-basin hopping method:
        # res = minimize(
        #     cost, [1, 1, 1, 1], args=(i, normalized_counts), method="Nelder-Mead"
        # )

        res = sp.optimize.basinhopping(
            state_cost,
            x0=[1, 1, 1, 1],
            minimizer_kwargs={
                "args": (i, normalized_counts, measurement_densities),
                #'options':{"disp":True}
            },
            niter=BASIN_HOPPING_N_ITER_STATE,
        )
        tomography_results.append(test_density(res.x))

    tomo_results = np.array(tomography_results)
    tomo_results = normalize_and_chop(tomo_results)
    purities = np.array([np.trace(density @ density) for density in tomo_results])

    return (tomo_results, purities, res.fun)


def get_constraint_from_chi(chi, basis=PAULIS):
    """
    Calculates the trace-preserving constraint for the chi matrix.

    Parameters:
        chi (numpy.ndarray): Chi matrix.
        basis (list): List of process basis matrices.
    Returns:
        numpy.ndarray: Matrix that is supposed to be the identity.
    """
    sum_ = np.zeros((2, 2), dtype=np.complex128)
    for i in range(4):
        for j in range(4):
            sum_ += chi[i, j].conj() * basis[i] @ basis[j]
    return sum_


def test_chi(params):
    """
    Calculates a chi matrix from a set of real parameters. This follows a lower-triangular parameterization.

    Parameters:
        chi (numpy.ndarray): Chi matrix.
        basis (list): List of process basis matrices.
    Returns:
        numpy.ndarray: Matrix that is supposed to be the identity.
    """
    L = np.array(
        [
            [params[0], 0, 0, 0],
            [params[4] + 1.0j * params[5], params[1], 0, 0],
            [params[10] + 1.0j * params[11], params[6] + 1.0j * params[7], params[2], 0],
            [
                params[14] + 1.0j * params[15],
                params[12] + 1.0j * params[13],
                params[8] + 1.0j * params[9],
                params[3],
            ],
        ],
        dtype=np.complex128,
    )

    return L @ (L.conj().T)


def normalize_and_chop_chi(chi):
    """
    Normalizes and chops the chi matrix to remove numerical artifacts.
    This is done by setting small values to zero and normalizing the trace.

    Parameters:
        chi (numpy.ndarray): Chi matrix.
    Returns:
        numpy.ndarray: Normalized and chopped chi matrix.
    """
    chi /= np.trace(chi)
    real_part = np.real(chi)
    imag_part = np.imag(chi)
    chopped_real = np.where(np.abs(real_part) < 1e-8, 0, real_part)
    chopped_imag = np.where(np.abs(imag_part) < 1e-8, 0, imag_part)
    ret = chopped_real + 1.0j * chopped_imag
    return ret


def get_process_output(chi, rho, basis=PAULIS):
    """
    Calculates the expected output state of a channel E[rho]
    for a given input state rho.

    Parameters:
        chi (numpy.ndarray): Chi matrix characterizing the channel E.
        rho (numpy.ndarray): Input density matrix.
    Returns:
        numpy.ndarray: Channel output density matrix.
    """
    expected_density = np.sum(
        [chi[i, j] * basis[i] @ rho @ basis[j] for i in range(4) for j in range(4)],
        axis=0,
        dtype=np.complex128,
    )
    expected_density /= np.trace(expected_density)
    return expected_density


def chi_cost(params, input_states, measurements, meas_counts):
    """
    Calculates chi cost function for process tomography.
    Compares the expected counts from the chi matrix with the measured counts.
    The cost function is the log likelihood of a Poissonian distribution.
    The cost function is calculated as (expected - actual)^2/(2*expected) for each measurement and sums them up.

    Parameters:
        params (numpy.ndarray): Real valued array characterizing the chi matrix.
        input_states (numpy.ndarray): Input density matrices.
        measurements (numpy.ndarray): Measurement density matrices for the measurements you performed.
        meas_counts (numpy.ndarray): Measured counts.
    Returns:
        float: Cost value.
    """
    # TODO: Make this not global
    global pre_computed_Ei_rho_Ej

    chi = test_chi(params)

    expected_densities = [
        np.dot(chi.flatten(), pre_computed_Ei_rho_Ej[i].reshape(16, 4)).reshape(2, 2)
        for i in range(len(input_states))
    ]
    # expected_densities = np.array([get_process_output(chi,state) for state in input_states],dtype=np.complex128)

    # Multiply expected densities matrix of shape (6,2,2) (i.e, xjk) with the input density matrices (6,2,2) (i.e, (ykl)) pairwise
    # to get an output matrix of (6,6,2,2) (i.e, xyjl), then trace over the resulting densities (i.e, axes 2 and 3, which are j and l)
    expected_counts = np.trace(
        np.einsum("xjk,ykl->xyjl", measurements, expected_densities),
        axis1=2,
        axis2=3,
        dtype=np.float128,
    )

    # # The line above basically performs this commented out section
    # expected_counts = np.zeros((6,6),dtype=np.complex128)
    # for i in range(6):
    #     for j in range(6):
    #         expected_counts[i,j] = np.trace(expected_densities[i]@measurements[j],dtype=float)

    # Clip The values so the cost function doesn't go crazy
    expected_counts = np.clip(expected_counts, 1e-10, None)

    diff_counts = expected_counts.flatten() - meas_counts.flatten()

    # cost = np.abs(np.sum((expected_counts-meas_counts.flatten())**2/(2*expected_counts)))

    # Calculate a matrix product (X.T@Y@X). The diagonal matrix adds the normalizing factor 1/2*expected_counts,
    # and the second X makes it squared
    # Equivalent to the line above
    cost = np.abs(diff_counts.T @ np.diag(1 / (2 * expected_counts.flatten())) @ diff_counts)

    return cost


def compute_Ei_rho_Ej_mat(input_densities, basis):
    """
    Used in calculation of chi_cost function.

    We can precompute output density matrices (before the sum) since we know the input states
    and the paulis are fixed, then we can apply the chi matrix to this later
    This gives a matrix like:
    [[E_1@rho@E_1, E_1@rho@E_2, E_1@rho@E_3, E_1@rho@E_4]
    [E_2@rho@E_1, E_2@rho@E_2, E_2@rho@E_3, E_2@rho@E_4]
    [E_3@rho@E_1, E_3@rho@E_2, E_3@rho@E_3, E_3@rho@E_4]
    [E_4@rho@E_1, E_4@rho@E_2, E_4@rho@E_3, E_4@rho@E_4]]

    Parameters:
        input_densities (numpy.ndarray): Input density matrices.
        basis (list): List of process basis matrices.
    Returns:
        numpy.ndarray: Precomputed matrix.
    """

    # TODO: Make this not global
    global pre_computed_Ei_rho_Ej

    pre_computed_Ei_rho_Ej = np.zeros(
        (len(input_densities), len(basis), len(basis), 2**NUM_QUBITS, 2**NUM_QUBITS),
        dtype=np.complex128,
    )
    for i in range(len(input_densities)):
        for m in range(len(basis)):
            for n in range(len(basis)):
                pre_computed_Ei_rho_Ej[i, m, n, 0:2, 0:2] = basis[m] @ input_densities[i] @ basis[n]

    return pre_computed_Ei_rho_Ej


def trace_preserving_constraint(params, basis=PAULIS):
    """
    Calculates the trace-preserving constraint for the chi matrix.

    Parameters:
        params (numpy.ndarray): Real-valued array parameterizing the chi matrix.
        basis (list): List of process basis matrices.
    Returns:
        numpy.ndarray: Trace-preserving constraint values (supposed to be equal to 0).
    """
    chi = test_chi(params)

    # The trace preserving constraint swaps the order of the basis matrices here,
    # so we use the conjugate transpose of the chi matrix (since chi is Hermitian,
    # we just take the transpose)

    S = sum(chi[j, i] * basis[i] @ basis[j] for i in range(4) for j in range(4))
    return [S[0, 0] - 1, S[1, 1] - 1, S[0, 1].imag, S[1, 0].real]


def get_rand_start_chi_parameters():
    """
    Get a random set of real-valued parameters for the chi matrix.

    Returns:
        numpy.ndarray: Real-valued parameters for chi matrix.
    """
    # TODO: Make this for generic number of qubits.
    rand_real = np.random.rand(16)
    rand_complex = np.random.rand(16)
    rand_start = rand_real + 1.0j * rand_complex
    return rand_start


def process_tomography(input_counts, output_counts, measurement_densities, basis=PAULIS):
    """
    Process tomography from measured counts. Input counts are used to calculate input states, and output counts are compared directly with expected process output counts.

    Parameters:
        input_counts (numpy.ndarray): Measured counts for your input states.
        output_counts (numpy.ndarray): Measured counts for the process output.
        measurement_densities (numpy.ndarray): Measurement density matrices for the measurements you performed.
        basis (list): List of process basis matrices.
    Returns:
        np.ndarray: Chi matrix characterizing the process.
        float: Final cost function value.
    """

    # TODO: Make this not global
    global pre_computed_Ei_rho_Ej

    normalized_output_counts = normalize_counts(output_counts)

    input_densities, _,_ = state_tomography(input_counts, measurement_densities)

    compute_Ei_rho_Ej_mat(input_densities, basis)

    ## If we want a random start:
    # rand_start = get_rand_start_chi_parameters()

    # Make random start a complex identity matrix, flattened and separated out the real and complex part (view(np.float64))
    rand_start = np.eye(4, dtype=complex).flatten().view(np.float64)
    constraints = [
        {"type": "eq", "fun": trace_preserving_constraint},
    ]

    # Debugging for constraints
    # The rank of the jacobian of your constraint function tells you how many constraints you actually need.
    # In our case, we only need 4 constraints instead of 8 (4 complex matrix elements),
    # since the trace preserving constraint is Hermitian.
    # So we can make sure that the diagonal elements are both 1, and we can choose either the
    # imaginary/real part of the off-diagonal elements to be 0.
    # from scipy.optimize import approx_fprime

    # jac = approx_fprime(rand_start, trace_preserving_constraint, 1e-8)
    # print("Constraint Jacobian shape:", jac.shape)
    # print("Rank of Jacobian:", np.linalg.matrix_rank(jac))

    # Different minimization methods
    # res = minimize(chi_cost, np.ones(16)*0.1, args=(input_densities, normalized_output_counts), constraints=constraints, method='COBYLA', options={"maxiter":10000})
    # res = minimize(chi_cost, x0=rand_start, args=(input_tomo_results,input_densities, normalized_output_counts), constraints=constraints, options={"maxiter":100000,"verbose":2})
    res = sp.optimize.basinhopping(
        chi_cost,
        x0=rand_start,
        minimizer_kwargs={
            "args": (input_densities, measurement_densities, normalized_output_counts),
            "constraints": constraints,
            #"options": {"disp": True},
        },
        niter=BASIN_HOPPING_N_ITER_PROCESS,
    )
    output_chi = test_chi(res.x)

    # Normalize and remove small values
    final_chi = normalize_and_chop([output_chi])[0]
    return final_chi, res.fun


if __name__ == "__main__":
    parser = argparse.ArgumentParser("Process tomography script for my data :^)")
    parser.add_argument("-f", "--file", help="File to load data from")
    parser.add_argument("--num_qubits", type=int, default=1, help="Number of qubits")
    parser.add_argument(
        "--basin_hopping_n_iter_state",
        type=int,
        default=10,
        help="Number of iterations for basin hopping in state tomography",
    )
    parser.add_argument(
        "--basin_hopping_n_iter_process",
        type=int,
        default=20,
        help="Number of iterations for basin hopping in process tomography",
    )
    args = parser.parse_args()

    # Setting globals TODO: Make this not use globals?
    BASIN_HOPPING_N_ITER_PROCESS = args.basin_hopping_n_iter_process
    BASIN_HOPPING_N_ITER_STATE = args.basin_hopping_n_iter_state
    NUM_QUBITS = args.num_qubits

    # Data importing stuff:
    with open(args.file, "r") as f:
        data = json.load(f)

    keys = list(data.keys())[5:]

    # Calculating background ratio from my data.
    # background_ratio = data["background_ratio"]
    # background_counts = data["background_counts"]
    # background_ratio = np.float128(background_counts[0][1] / data["shot_rate"])

    # This is weird funky stuff because I messed up the measurements, ignore
    counts = []
    for key in keys:
        if key == "RL":
            data_arr = np.array(data["RR"])
        if key == "LR":
            data_arr = np.array(data["LL"])
        if key == "LL":
            data_arr = np.array(data["LR"])
        if key == "RR":
            data_arr = np.array(data["RL"])
        else:
            data_arr = np.array(data[key])

        counts.append(np.mean(data_arr[:, 1]))

    input_counts_H = [1, 0, 0.5, 0.5, 0.5, 0.5]
    input_counts_V = [0, 1, 0.5, 0.5, 0.5, 0.5]
    input_counts_D = [0.5, 0.5, 1, 0, 0.5, 0.5]
    input_counts_A = [0.5, 0.5, 0, 1, 0.5, 0.5]
    input_counts_R = [0.5, 0.5, 0.5, 0.5, 1, 0]
    input_counts_L = [0.5, 0.5, 0.5, 0.5, 0, 1]

    input_counts = np.array(
        [
            input_counts_H,
            input_counts_V,
            input_counts_D,
            input_counts_A,
            input_counts_R,
            input_counts_L,
        ]
    ).reshape(6, 6)

    output_counts = np.array(counts).reshape(6, 6)

    H = np.array([1, 0])
    V = np.array([0, 1])
    D = 1 / np.sqrt(2) * np.array([1, 1])
    A = 1 / np.sqrt(2) * np.array([1, -1])
    R = 1 / np.sqrt(2) * np.array([1, 1.0j])
    L = 1 / np.sqrt(2) * np.array([1, -1.0j])

    input_states = [H, V, D, A, R, L]
    densities = []
    for state in input_states:
        density = np.outer(state, np.conj(state))
        densities.append(density)
    measurement_densities = np.array(densities)

    output_chi, cost = process_tomography(input_counts, output_counts, measurement_densities)

    # Chi matrix for the identity channel
    identity_channel = np.zeros((4,4), dtype=np.complex128)
    identity_channel[0,0] = 1

    print(sym.pprint(sym.Matrix(output_chi), use_unicode=False))
    print("Cost function value: ", cost)
    print("Fidelity: ", calculate_fidelity(output_chi, identity_channel))

    ## For saving chi matrix to file
    # import pickle
    # from pathlib import Path
    # import datetime
    # 
    # data_path = Path("./process_tomo_results/")
    # now = datetime.datetime.now()
    # data_path = (
    #     data_path
    #     / f"{now.strftime('%Y-%m-%dT%H-%M-%SZ')}_{args.file.strip(".json").strip("/")}_scipy_process_tomo_results.pkl"
    # )
    # Path(data_path).mkdir(exist_ok=True)
    # with open(Path(data_path), "wb") as f:
    #     pickle.dump(output_chi, f)
