from itertools import product

import numpy as np
import pytest
from scipy.stats import unitary_group

import QuantumTomography as qLib

# Mirrors the "measurement_states" convention used in ExampleFiles/*.json:
# a single dict, shared across all qubits, mapping a basis label to the
# single-qubit state it projects onto.
STANDARD_BASIS = {
    "H": np.array([1, 0], dtype=complex),
    "V": np.array([0, 1], dtype=complex),
    "D": np.array([1, 1], dtype=complex) / np.sqrt(2),
    "A": np.array([1, -1], dtype=complex) / np.sqrt(2),
    "R": np.array([1, 1j], dtype=complex) / np.sqrt(2),
    "L": np.array([1, -1j], dtype=complex) / np.sqrt(2),
}
BASIS_NAMES = list(STANDARD_BASIS.keys())

# H/V, D/A and R/L are each an orthonormal basis, and the three bases are
# mutually unbiased (any state from one basis has overlap probability 1/2
# with either state of the other two bases).
MUB_GROUPS = [["H", "V"], ["D", "A"], ["R", "L"]]

N_COUNTS = 3000
FIDELITY_THRESHOLD = 0.95
N_TRIALS = 5


def random_mub_basis():
    """Return a randomly rotated version of H/V, D/A, R/L that is still mutually unbiased.

    A global unitary rotation preserves inner products, so applying the same
    random unitary to all 6 standard states keeps the 3 rotated pairs
    mutually unbiased. This is the plain quantum-mechanical state; it is not
    yet in the form the tomography library expects (see to_measurement_state).
    """
    u = unitary_group.rvs(2)
    return {label: u @ state for label, state in STANDARD_BASIS.items()}


def assert_mutually_unbiased(basis, atol=1e-6):
    for label_a, label_b in MUB_GROUPS:
        assert abs(np.vdot(basis[label_a], basis[label_b])) < atol

    for i in range(len(MUB_GROUPS)):
        for j in range(i + 1, len(MUB_GROUPS)):
            for label_a in MUB_GROUPS[i]:
                for label_b in MUB_GROUPS[j]:
                    overlap_prob = abs(np.vdot(basis[label_a], basis[label_b])) ** 2
                    assert abs(overlap_prob - 0.5) < atol


def to_measurement_state(state):
    """Convert a pure state to the form the tomography measurement matrix expects.

    Global phase is physically unobservable, but the library's internal
    projector construction (TomoClass.filter_data) implicitly assumes each
    measurement's first amplitude is real. TODO: We should fix this. Standard
    polarization states (H/V/D/A/R/L) already satisfy this; an arbitrarily
    rotated state may not, so its global phase is fixed here before it is
    placed in the measurement matrix.
    """
    return qLib.removeGlobalPhase(state)


def build_measurements(measurement_states, basis_names, num_qubits):
    """Build measurement settings the same way ExampleFiles/*.json data entries do.

    Each measurement setting is a "basis" list of `num_qubits` labels (one per
    qubit, e.g. ["D", "R"]), looked up in a single measurement_states dict
    shared across all qubits (see e.g. ExampleFiles/bell_state_example.json).
    Every combination of basis_names across the qubits is covered. The result
    has shape (len(basis_names)**num_qubits, num_qubits, 2): one row per
    measurement setting, one (alpha, beta) pair per qubit within that
    setting. This is the plain physical representation, not yet converted to
    the flat layout StateTomography expects (see to_tomo_input_measurements).
    """
    bases = list(product(basis_names, repeat=num_qubits))
    measurements = np.zeros((len(bases), num_qubits, 2), dtype=complex)
    for row, basis in enumerate(bases):
        for qubit_index, label in enumerate(basis):
            measurements[row, qubit_index] = measurement_states[label]
    return measurements


def to_tomo_input_measurements(measurements):
    """Convert measurement states into the format that the library expects.

    Flatten (M, num_qubits, 2) measurement settings into the (M, 2*num_qubits)
    alpha/beta-per-qubit layout StateTomography expects, fixing each state's
    global phase along the way (see to_measurement_state).
    """
    num_measurements, num_qubits, _ = measurements.shape
    flattened = np.empty((num_measurements, 2 * num_qubits), dtype=complex)
    for i in range(num_measurements):
        for q in range(num_qubits):
            flattened[i, 2 * q : 2 * q + 2] = to_measurement_state(measurements[i, q])
    return flattened


@pytest.fixture
def rng(request):
    # Ties our Generator to pytest-randomly's --randomly-seed, so a failing
    # run (including the simulated shot noise) can be reproduced exactly via
    # `pytest --randomly-seed=<seed>`.
    seed = request.config.getoption("randomly_seed")
    return np.random.default_rng(seed)


def simulate_counts(measurements, rho, rng, n_counts=N_COUNTS):
    # TODO: Move this to utilities?
    # Single-detector-per-qubit setup: one Poisson-sampled photon count per measurement setting.
    num_measurements, num_qubits, _ = measurements.shape
    counts = np.zeros(num_measurements)
    for i in range(num_measurements):
        state = np.array([1], dtype=complex)
        for q in range(num_qubits):
            state = np.kron(state, measurements[i, q])
        prob = np.real(np.vdot(state, rho @ state))
        counts[i] = rng.poisson(n_counts * np.clip(prob, 0, 1))
    return counts


@pytest.mark.parametrize("num_qubits", [1, 2])
@pytest.mark.parametrize("_trial", range(N_TRIALS))
def test_state_tomography_standard_basis(num_qubits, _trial, rng):
    rho = qLib.random_density_state(num_qubits)

    measurements = build_measurements(STANDARD_BASIS, BASIS_NAMES, num_qubits)
    counts = simulate_counts(measurements, rho, rng)

    tomo_obj = qLib.Tomography()
    [rho_approx, _, _] = tomo_obj.StateTomography(to_tomo_input_measurements(measurements), counts)

    assert qLib.fidelity(rho, rho_approx) > FIDELITY_THRESHOLD


@pytest.mark.parametrize("num_qubits", [1, 2])
@pytest.mark.parametrize("_trial", range(N_TRIALS))
def test_state_tomography_random_mub(num_qubits, _trial, rng):
    rho = qLib.random_density_state(num_qubits)

    basis = random_mub_basis()
    assert_mutually_unbiased(basis)

    measurements = build_measurements(basis, BASIS_NAMES, num_qubits)
    counts = simulate_counts(measurements, rho, rng)

    tomo_obj = qLib.Tomography()
    [rho_approx, _, _] = tomo_obj.StateTomography(to_tomo_input_measurements(measurements), counts)

    assert qLib.fidelity(rho, rho_approx) > FIDELITY_THRESHOLD
