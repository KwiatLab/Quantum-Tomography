from functools import reduce
from itertools import product

import numpy as np
import scipy
from hypothesis import assume
from hypothesis import example
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays

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
TOLERANCE = 1e-9  # seems to work out to 1e-15


@st.composite
def norm(draw, loc=0, scale=1):
    seed = draw(st.integers(0, 2**31 - 1))
    return scipy.stats.norm.rvs(loc=loc, scale=scale, random_state=seed)


@st.composite
def poisson(draw, mu=1, size=1):
    seed = draw(st.integers(0, 2**31 - 1))
    return scipy.stats.poisson.rvs(mu=mu, size=size, random_state=seed)


@st.composite
def binormal_complex(draw):
    a = draw(norm())
    b = draw(norm())
    return a + b * 1j


@st.composite
def ginibre_matrix(draw, dim=2):
    return draw(arrays(dtype=complex, shape=((dim, dim)), elements=binormal_complex(), fill=st.nothing()))


@st.composite
def density_matrix(draw, n_qubits=1):
    _gm = draw(ginibre_matrix(2**n_qubits))
    density = np.matmul(_gm, _gm.T.conj())

    trace = np.trace(density)
    assume(trace >= TOLERANCE)

    return density / trace


@st.composite
def unitary_matrix(draw, dim=2):
    seed = draw(st.integers(0, 2**31 - 1))
    return scipy.stats.unitary_group.rvs(dim, random_state=seed)


@example(STANDARD_BASIS)
@st.composite
def mutually_unbiased_basis(draw):
    """Return a randomly rotated version of H/V, D/A, R/L that is still mutually unbiased.

    A global unitary rotation preserves inner products, so applying the same
    random unitary to all 6 standard states keeps the 3 rotated pairs
    mutually unbiased. This is the plain quantum-mechanical state; it is not
    yet in the form the tomography library expects (see to_measurement_state).
    """
    u = draw(unitary_matrix())
    basis = {label: u @ state for label, state in STANDARD_BASIS.items()}

    # TODO: Vectorize for a single `assume()` call
    for label_a, label_b in MUB_GROUPS:
        assume(abs(np.vdot(basis[label_a], basis[label_b])) < TOLERANCE)

    for i in range(len(MUB_GROUPS)):
        for j in range(i + 1, len(MUB_GROUPS)):
            for label_a in MUB_GROUPS[i]:
                for label_b in MUB_GROUPS[j]:
                    overlap_prob = abs(np.vdot(basis[label_a], basis[label_b])) ** 2
                    assume(abs(overlap_prob - 0.5) < TOLERANCE)

    return basis


def build_measurements(measurement_states, num_qubits):
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
    return np.array(
        [
            [measurement_states[label] for label in measurement_setting]
            for measurement_setting in product(BASIS_NAMES, repeat=num_qubits)
        ],
        dtype=complex,
    )


def remove_global_phase(states):
    """Convert a pure state to the form the tomography measurement matrix expects.

    Global phase is physically unobservable, but the library's internal
    projector construction (TomoClass.filter_data) implicitly assumes each
    measurement's first amplitude is real. TODO: We should fix this. Standard
    polarization states (H/V/D/A/R/L) already satisfy this; an arbitrarily
    rotated state may not, so its global phase is fixed here before it is
    placed in the measurement matrix.
    """
    norms = np.sqrt(np.sum(states.conj() * states, axis=1, keepdims=True))
    states_norm = states / norms
    phases = np.angle(states_norm[:, 0:1])  # (N, 1)
    phase_factors = np.exp(1j * phases)
    return states_norm / phase_factors


def to_tomo_input_measurements(measurements):
    """Convert measurement states into the format that the library expects.

    Flatten (M, num_qubits, 2) measurement settings into the (M, 2*num_qubits)
    alpha/beta-per-qubit layout StateTomography expects, fixing each state's
    global phase along the way (see remove_global_phase).
    """
    num_measurements, num_qubits, _ = measurements.shape
    reshaped = measurements.reshape(-1, 2)
    processed = remove_global_phase(reshaped)
    return processed.reshape(num_measurements, 2 * num_qubits)


@st.composite
def simulate_counts(draw, rho, measurements, n_counts=N_COUNTS):
    # Single-detector-per-qubit setup: one Poisson-sampled photon count per measurement setting.
    states = np.array([reduce(np.kron, measurement_setting) for measurement_setting in measurements], dtype=complex)

    # Compute conj(state) @ (rho @ state) for each measurement setting in a single operation.
    probs = np.einsum("mi,ij,mj->m", states.conj(), rho, states).real.clip(0, 1)

    counts = draw(poisson(mu=n_counts, size=probs.size))
    return probs * counts


@st.composite
def state_tomo_data(draw):
    n_qubits = draw(st.integers(1, 2))
    rho = draw(density_matrix(n_qubits))

    basis = draw(mutually_unbiased_basis())
    measurements = build_measurements(basis, n_qubits)

    counts = draw(simulate_counts(rho, measurements))
    return rho, measurements, counts


@given(data=state_tomo_data())
@settings(max_examples=N_TRIALS)
def test_state_tomography(data):
    rho, measurements, counts = data

    tomo_obj = qLib.Tomography()
    [rho_approx, _, _] = tomo_obj.StateTomography(to_tomo_input_measurements(measurements), counts)

    assert qLib.fidelity(rho, rho_approx) > FIDELITY_THRESHOLD
