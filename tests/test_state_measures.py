"""First-principles tests of the state measures and representations used by tomography.

Expected values come from closed-form results, never from the library itself.
Some of these fail on purpose: they document open bugs listed in TODO.md.
"""

import numpy as np
import numpy.testing as npt
import scipy
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from test_state_tomography import density_matrix
from test_state_tomography import ginibre_matrix
from test_state_tomography import unitary_matrix

import QuantumTomography as qLib

N_EXAMPLES = 50
FULL_RANK = 1e-6  # smallest eigenvalue of a state that counts as full rank
PSI_MINUS = np.array([0, 1, -1, 0], dtype=complex) / np.sqrt(2)


def projector(state):
    return np.outer(state, state.conj())


def werner(p):
    """Return p |psi-><psi-| + (1 - p) I/4, which is entangled iff p > 1/3."""
    return p * projector(PSI_MINUS) + (1 - p) * np.eye(4) / 4


@st.composite
def pure_state(draw, dim=4):
    """Draw a Haar-random state vector."""
    return draw(unitary_matrix(dim))[:, 0]


@st.composite
def n_qubits_and_states(draw, n_states=1):
    """Draw full-rank states of one or two qubits.

    Hypothesis sometimes repeats a seed, which makes density_matrix() return a
    rank-deficient state; density2t crashes on those (see TODO.md), which
    tests/test_tomography_physics.py covers with pure states on purpose.
    """
    n_qubits = draw(st.integers(1, 2))
    states = [draw(density_matrix(n_qubits)) for _ in range(n_states)]
    for rho in states:
        assume(np.linalg.eigvalsh(rho)[0] > FULL_RANK)
    return n_qubits, states


@given(data=n_qubits_and_states(n_states=2))
@settings(max_examples=N_EXAMPLES)
def test_fidelity_matches_uhlmann_formula(data):
    """F = (Tr sqrt(sqrt(rho) sigma sqrt(rho)))^2, and it is symmetric."""
    _, (rho, sigma) = data
    sqrt_rho = scipy.linalg.sqrtm(rho)
    expected = np.sum(np.sqrt(np.linalg.eigvals(sqrt_rho @ sigma @ sqrt_rho).real.clip(0))) ** 2

    npt.assert_allclose(qLib.fidelity(rho, sigma), expected, atol=1e-6)
    npt.assert_allclose(qLib.fidelity(sigma, rho), expected, atol=1e-6)


@given(a=pure_state(), b=pure_state())
@settings(max_examples=N_EXAMPLES)
def test_fidelity_of_pure_states_is_squared_overlap(a, b):
    expected = abs(np.vdot(a, b)) ** 2

    npt.assert_allclose(qLib.fidelity(a, b), expected, atol=1e-9)
    npt.assert_allclose(qLib.fidelity(projector(a), b), expected, atol=1e-9)
    npt.assert_allclose(qLib.fidelity(projector(a), projector(b)), expected, atol=1e-6)


@given(u=unitary_matrix(4), weights=st.lists(st.floats(0.01, 1), min_size=4, max_size=4))
@settings(max_examples=N_EXAMPLES)
def test_purity_and_entropy_depend_only_on_spectrum(u, weights):
    spectrum = np.array(weights) / np.sum(weights)
    rho = u @ np.diag(spectrum) @ u.conj().T

    npt.assert_allclose(qLib.purity(rho), np.sum(spectrum**2), atol=1e-9)
    npt.assert_allclose(qLib.linear_entropy(rho), 1 - np.sum(spectrum**2), atol=1e-9)
    npt.assert_allclose(qLib.entropy(rho), -np.sum(spectrum * np.log2(spectrum)), atol=1e-9)


@given(psi=pure_state())
@settings(max_examples=N_EXAMPLES)
def test_entanglement_of_pure_states(psi):
    """For a|00> + b|01> + c|10> + d|11>, concurrence and negativity are both 2 |ad - bc|."""
    a, b, c, d = psi
    expected = 2 * abs(a * d - b * c)
    rho = projector(psi)

    npt.assert_allclose(qLib.concurrence(rho), expected, atol=1e-6)
    npt.assert_allclose(qLib.tangle(rho), expected**2, atol=1e-6)
    npt.assert_allclose(qLib.negativity(rho), expected, atol=1e-6)


@given(p=st.floats(0, 1))
@settings(max_examples=N_EXAMPLES)
def test_entanglement_of_werner_states(p):
    """Concurrence and negativity are both max(0, (3p - 1) / 2); the state is separable for p <= 1/3."""
    expected = max(0.0, (3 * p - 1) / 2)

    npt.assert_allclose(qLib.concurrence(werner(p)), expected, atol=1e-6)
    npt.assert_allclose(qLib.negativity(werner(p)), expected, atol=1e-6)


@given(data=n_qubits_and_states())
@settings(max_examples=N_EXAMPLES)
def test_t_parametrization_round_trip(data):
    """MLE optimizes over t values, so density -> t -> density must be the identity."""
    _, (rho,) = data

    npt.assert_allclose(qLib.t_to_density(qLib.density2t(rho)), rho, atol=1e-7)


@given(p=st.floats(0, 1), noise=ginibre_matrix(4))
@settings(max_examples=N_EXAMPLES)
def test_make_positive_keeps_physical_states(p, noise):
    """Rounding noise on a state with repeated eigenvalues (any Werner state) must not change it."""
    rho = werner(p) + 1e-15 * (noise + noise.conj().T)

    npt.assert_allclose(qLib.make_positive(rho), rho, atol=1e-9)


@given(data=n_qubits_and_states())
@settings(max_examples=N_EXAMPLES)
def test_state_is_rebuilt_from_stokes_parameters(data):
    """Linear tomography relies on rho = sum_i S_i sigma_i / 2^n with S_i = Tr(rho sigma_i)."""
    n_qubits, (rho,) = data
    stokes = qLib.get_stokes_parameters(rho)
    paulis = qLib.generalized_pauli_basis(n_qubits)

    rebuilt = np.tensordot(stokes, paulis, axes=1) / 2**n_qubits

    npt.assert_allclose(rebuilt, rho, atol=1e-9)


@given(psi=pure_state(dim=2), qwp=st.floats(0, np.pi), hwp=st.floats(0, np.pi), flip=st.booleans())
@settings(max_examples=N_EXAMPLES)
def test_wave_plate_basis_is_projection_seen_at_pbs(psi, qwp, hwp, flip):
    """Light passes the QWP, then the HWP, then a PBS that transmits H (or V if flipped)."""
    after_plates = qLib.halfWavePlate(hwp) @ qLib.quarterWavePlate(qwp) @ psi
    p_transmit = abs(after_plates[1 if flip else 0]) ** 2

    basis = qLib.getWavePlateBasis(qwp, hwp, flipPBS=flip)

    npt.assert_allclose(abs(np.vdot(basis[0], psi)) ** 2, p_transmit, atol=1e-9)
    npt.assert_allclose(abs(np.vdot(basis[1], psi)) ** 2, 1 - p_transmit, atol=1e-9)
