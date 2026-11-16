"""First-principles tests of state tomography and its experimental corrections.

Each test prepares a known state, computes the counts an ideal experiment would
record from the Born rule (without using the library), applies the imperfection
under test, and checks that the library recovers the state that was prepared.

Several of these fail on purpose: they document open bugs listed in TODO.md.
"""

import json
import tempfile
from functools import reduce
from itertools import product
from pathlib import Path

import numpy as np
import pytest
from hypothesis import Phase
from hypothesis import assume
from hypothesis import given
from hypothesis import settings
from hypothesis import strategies as st
from hypothesis.extra.numpy import arrays
from test_state_tomography import BASIS_NAMES
from test_state_tomography import MUB_GROUPS
from test_state_tomography import N_TRIALS
from test_state_tomography import STANDARD_BASIS
from test_state_tomography import density_matrix
from test_state_tomography import unitary_matrix

import QuantumTomography as qLib

RATE = 1000.0  # coincidences per measurement for a projector with probability 1
JSON_RATE = 1e5  # JSON counts are whole numbers; high enough that rounding is negligible
EXACT_FIDELITY = 0.999  # noise-free data
FULL_RANK = 1e-6  # smallest eigenvalue of a state that counts as full rank
GOOD_FIT = 1.0  # fval (chi-squared) of data the model can describe exactly
HEDGE_BETA = 1e-3
N_SIGMA = 5

# Shrinking is switched off: every example is a full fit (up to seconds when the fit goes wrong), and the
# states come from opaque seeds, so shrinking a failure costs minutes and does not make the example simpler.
tomography_settings = settings(max_examples=N_TRIALS, phases=[Phase.explicit, Phase.reuse, Phase.generate])

PHI_PLUS = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
PHI_MINUS = np.array([1, 0, 0, -1], dtype=complex) / np.sqrt(2)
PSI_PLUS = np.array([0, 1, 1, 0], dtype=complex) / np.sqrt(2)
PSI_MINUS = np.array([0, 1, -1, 0], dtype=complex) / np.sqrt(2)
BELL_STATES = {"phi+": PHI_PLUS, "phi-": PHI_MINUS, "psi+": PSI_PLUS, "psi-": PSI_MINUS}


def projector(state):
    return np.outer(state, state.conj()) / np.vdot(state, state).real


def werner(p):
    """Return p |psi-><psi-| + (1 - p) I/4, which is entangled iff p > 1/3."""
    return p * projector(PSI_MINUS) + (1 - p) * np.eye(4) / 4


def standard_measurements(n_qubits, n_det):
    """Build the measurement matrix StateTomography expects: one (alpha, beta) pair per qubit, flattened.

    With one detector per qubit every setting is a single projection, so all six
    states are needed. With two detectors the orthogonal state is recorded at
    the same time, so one state from each of H/V, D/A and R/L is enough.
    """
    labels = BASIS_NAMES if n_det == 1 else [group[0] for group in MUB_GROUPS]
    return np.array(
        [np.concatenate([STANDARD_BASIS[label] for label in setting]) for setting in product(labels, repeat=n_qubits)]
    )


def detected_states(setting, n_det):
    """Return every product state detected for one measurement setting.

    With two detectors per qubit each qubit is split into the given state and
    the one orthogonal to it. Outcomes are ordered as binary numbers with qubit
    0 as the most significant bit (for two qubits: 11, 12, 21, 22).
    """
    qubit_states = [state / np.linalg.norm(state) for state in setting.reshape(-1, 2)]
    if n_det == 1:
        return [reduce(np.kron, qubit_states)]
    pairs = [(state, np.array([state[1].conj(), -state[0].conj()])) for state in qubit_states]
    return [reduce(np.kron, choice) for choice in product(*pairs)]


def probabilities(rho, measurements, n_det):
    """Born rule <psi|rho|psi>, with shape (n_measurements, n_outcomes)."""
    return np.array(
        [[np.vdot(state, rho @ state).real for state in detected_states(setting, n_det)] for setting in measurements]
    )


def ideal_counts(rho, measurements, n_det, rate=RATE):
    counts = rate * probabilities(rho, measurements, n_det)
    return counts[:, 0] if n_det == 1 else counts


def tomography(rho, n_det=1, method="MLE", counts=None, conf=None, **kwargs):
    """Run the library on (by default ideal) counts from rho; returns (estimate, intensity, fval, tomo_obj)."""
    n_qubits = int(np.log2(rho.shape[0]))
    measurements = kwargs.pop("measurements", standard_measurements(n_qubits, n_det))
    if counts is None:
        counts = ideal_counts(rho, measurements, n_det)
    tomo_obj = qLib.Tomography(n_qubits)
    for key, value in (conf or {}).items():
        tomo_obj.conf[key] = value
    estimate, intensity, fval = tomo_obj.StateTomography(measurements, counts, method=method, **kwargs)
    return estimate, intensity, fval, tomo_obj


@st.composite
def state(draw, n_qubits=None):
    """Draw a random full-rank density matrix of one or two qubits.

    Hypothesis sometimes repeats a seed, which makes density_matrix() return a
    rank-deficient state. Those hit the pure-state bugs, which
    test_ideal_counts_recover_pure_state covers on purpose, so they are
    excluded here to keep every other test deterministic.
    """
    if n_qubits is None:
        n_qubits = draw(st.integers(1, 2))
    rho = draw(density_matrix(n_qubits))
    assume(np.linalg.eigvalsh(rho)[0] > FULL_RANK)
    return rho


@st.composite
def pure_state(draw):
    n_qubits = draw(st.integers(1, 2))
    return projector(draw(unitary_matrix(2**n_qubits))[:, 0])


def factors(low, high, size):
    return arrays(float, size, elements=st.floats(low, high))


# Reconstruction from ideal data


@pytest.mark.parametrize("method", ["MLE", "LINEAR"])
@pytest.mark.parametrize("n_det", [1, 2])
@given(rho=state())
@tomography_settings
def test_ideal_counts_recover_state(n_det, method, rho):
    estimate, intensity, _, _ = tomography(rho, n_det, method)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY
    # counts are rate * probability, so the fitted intensity is the rate whichever method is used
    assert intensity == pytest.approx(RATE, rel=1e-2)


@given(rho=pure_state())
@settings(tomography_settings, max_examples=5 * N_TRIALS)
def test_ideal_counts_recover_pure_state(rho):
    """Pure states lie on the boundary of state space; noise-free data from one must not crash the fit."""
    estimate, _, _, _ = tomography(rho)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY


@pytest.mark.parametrize("n_det", [1, 2])
def test_three_qubit_ghz_state(n_det):
    ghz = np.zeros(8, dtype=complex)
    ghz[[0, 7]] = 1 / np.sqrt(2)

    estimate, _, _, _ = tomography(projector(ghz), n_det)

    assert qLib.fidelity(ghz, estimate) > EXACT_FIDELITY


@given(rho=state(n_qubits=1), angle=st.floats(0, 2 * np.pi))
@tomography_settings
def test_global_phase_of_measurement_states_is_irrelevant(rho, angle):
    """|psi> and e^{i phi}|psi> are the same projector, so the same counts must give the same state."""
    measurements = standard_measurements(1, 1)
    counts = ideal_counts(rho, measurements, 1)

    estimate, _, _, _ = tomography(rho, counts=counts, measurements=measurements * np.exp(1j * angle))

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY


# Experimental corrections


@pytest.mark.parametrize("n_det", [1, 2])
@given(rho=state(n_qubits=2), data=st.data())
@tomography_settings
def test_drift_correction(n_det, rho, data):
    """If the pump power drifts, setting i records rate * intensity_i * probability."""
    measurements = standard_measurements(2, n_det)
    drift = data.draw(factors(0.5, 1.5, len(measurements)))
    counts = ideal_counts(rho, measurements, n_det) * (drift if n_det == 1 else drift[:, None])

    estimate, intensity, fval, _ = tomography(rho, n_det, counts=counts, intensities=drift)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY
    assert intensity == pytest.approx(RATE, rel=1e-2)
    assert fval < GOOD_FIT


@given(rho=state(), data=st.data())
@tomography_settings
def test_integration_time_normalization(rho, data):
    """Counting for twice as long doubles the counts; the estimate must not change."""
    n_qubits = int(np.log2(rho.shape[0]))
    measurements = standard_measurements(n_qubits, 1)
    times = data.draw(factors(1, 3, len(measurements)))
    counts = ideal_counts(rho, measurements, 1) * times

    estimate, _, _, _ = tomography(rho, counts=counts, time=times)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY


@given(rho=state(n_qubits=2), efficiency=factors(0.5, 1, 4))
@tomography_settings
def test_detector_efficiency_correction(rho, efficiency):
    """With two detectors per qubit, outcome k is recorded with relative efficiency efficiency[k]."""
    counts = ideal_counts(rho, standard_measurements(2, 2), 2) * efficiency

    estimate, _, fval, _ = tomography(rho, 2, counts=counts, efficiency=efficiency)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY
    assert fval < GOOD_FIT


@given(rho=state(n_qubits=1), leak=factors(0, 0.3, 2))
@tomography_settings
def test_crosstalk_correction(rho, leak):
    """A lossless beam splitter leaks a fraction of each port into the other.

    crosstalk[k, l] is the probability that a photon belonging in port l is
    detected in port k, so each column sums to one.
    """
    crosstalk = np.array([[1 - leak[0], leak[1]], [leak[0], 1 - leak[1]]])
    counts = RATE * probabilities(rho, standard_measurements(1, 2), 2) @ crosstalk.T

    estimate, _, _, _ = tomography(rho, 2, counts=counts, crosstalk=crosstalk)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY


@pytest.mark.parametrize("n_det", [1, 2])
@given(rho=state(n_qubits=2))
@tomography_settings
def test_accidental_correction(n_det, rho):
    """Uncorrelated singles S1, S2 give S1 * S2 * window / time accidental coincidences."""
    measurements = standard_measurements(2, n_det)
    singles_rate, window_ns, seconds = 2e5, 5.0, 2.0
    accidentals = singles_rate**2 * window_ns * 1e-9 / seconds  # 100 on top of at most 1000 real counts
    counts = ideal_counts(rho, measurements, n_det) + accidentals

    estimate, intensity, fval, _ = tomography(
        rho,
        n_det,
        counts=counts,
        conf={"DoAccidentalCorrection": 1},
        singles=np.full((len(measurements), 2 * n_det), singles_rate),
        time=np.full(len(measurements), seconds),
        window=np.full(1 if n_det == 1 else 4, window_ns),
    )

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY
    assert intensity == pytest.approx(RATE, rel=1e-2)
    assert fval < GOOD_FIT


# Estimators and derived quantities


@pytest.mark.parametrize("n_det", [1, 2])
def test_hedged_mle_with_small_beta_recovers_state(n_det):
    """With a negligible beta, hedged MLE is plain MLE.

    TODO: this uses one fixed state because the hedged fit converges poorly on
    random states (fidelity 0.94-1.0 on exact data, see TODO.md); once that is
    fixed, add "HMLE" to the methods of test_ideal_counts_recover_state.
    """
    rho = werner(0.6)

    estimate, intensity, _, _ = tomography(rho, n_det, "HMLE", conf={"Beta": HEDGE_BETA})

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY
    assert intensity == pytest.approx(RATE, rel=1e-2)


def test_hedged_mle_follows_add_beta_rule():
    """Hedging maximizes likelihood * det(rho)^beta (Blume-Kohout 2010).

    If N photons are all seen in H, the hedged estimate of the V population is
    beta / (N + 2 beta), where plain maximum likelihood gives exactly zero.
    """
    beta, n_photons = 0.5, 20.0
    h_state = projector(STANDARD_BASIS["H"])
    measurements = standard_measurements(1, 1)
    counts = ideal_counts(h_state, measurements, 1, rate=n_photons)

    estimate, _, _, _ = tomography(h_state, method="HMLE", counts=counts, conf={"Beta": beta})

    v_population = np.linalg.eigvalsh(estimate)[0]
    assert v_population == pytest.approx(beta / (n_photons + 2 * beta), rel=0.5)


def test_monte_carlo_error_bars_cover_true_values():
    p, rate, n_runs = 0.8, 5000.0, 12
    rho = werner(p)
    measurements = standard_measurements(2, 1)
    counts = np.random.default_rng(3).poisson(ideal_counts(rho, measurements, 1, rate)).astype(float)
    np.random.seed(3)  # noqa: NPY002 - the library resamples with the global numpy generator

    _, _, _, tomo_obj = tomography(rho, counts=counts)
    properties = {name: (mean, std) for name, mean, std in tomo_obj.getProperties(n_runs)}

    assert len(tomo_obj.mont_carlo_states) == n_runs + 1
    truth = {"intensity": rate, "concurrence": (3 * p - 1) / 2, "purity": (1 + 3 * p**2) / 4}
    for name, value in truth.items():
        mean, std = properties[name]
        assert std > 0
        assert abs(mean - value) < N_SIGMA * std + 0.01 * value, name


@pytest.mark.parametrize("name", BELL_STATES)
def test_bell_states_reach_tsirelson_bound(name):
    """Every maximally entangled state violates CHSH maximally: S = 2 sqrt(2)."""
    _, _, _, tomo_obj = tomography(projector(BELL_STATES[name]))

    settings_table = {row[0]: row[1] for row in tomo_obj.getBellSettings()}

    assert settings_table["s"] == pytest.approx(2 * np.sqrt(2), abs=1e-3)


# JSON data files


def json_data(rho, n_det):
    """Describe an ideal experiment in the format of ExampleFiles/*.json.

    counts is [singles for each qubit ..., lower-fold coincidences ...,
    n-fold coincidences]; only the last entry matters here. With two detectors,
    qubit q uses detector q for the first state of a pair and q + n_qubits for
    the second, and the outcomes of one setting are listed together.
    """
    n_qubits = int(np.log2(rho.shape[0]))
    if n_det == 1:
        settings_and_outcomes = [(setting, [0] * n_qubits) for setting in product(BASIS_NAMES, repeat=n_qubits)]
    else:
        settings_and_outcomes = [
            ([pair[o] for pair, o in zip(pairs, outcome, strict=True)], outcome)
            for pairs in product(MUB_GROUPS, repeat=n_qubits)
            for outcome in product([0, 1], repeat=n_qubits)
        ]

    data = []
    for labels, outcome in settings_and_outcomes:
        measured = reduce(np.kron, [STANDARD_BASIS[label] for label in labels])
        coincidences = round(JSON_RATE * np.vdot(measured, rho @ measured).real)
        datum = {"basis": list(labels), "integration_time": 1, "counts": [0] * (2**n_qubits - 2) + [coincidences]}
        if n_det > 1:
            datum["detectors_used"] = [q + o * n_qubits for q, o in enumerate(outcome)]
        data.append(datum)

    return {
        "n_qubits": n_qubits,
        "n_detectors_per_qubit": n_det,
        "n_measurements_per_qubit": len(BASIS_NAMES),
        "measurement_states": {"H": [1, 0], "V": [0, 1], "D": [1, 1], "A": [1, -1], "R": [1, "1j"], "L": [1, "-1j"]},
        "data": data,
        # TODO: a one-detector file without this key crashes on import (see TODO.md)
        **({"coincidence_window": [0]} if n_det == 1 else {}),
    }


def tomography_from_json(data):
    with tempfile.TemporaryDirectory() as directory:
        path = Path(directory) / "data.json"
        path.write_text(json.dumps(data))
        tomo_obj = qLib.Tomography()
        tomo_obj.import_data(str(path))
        return tomo_obj.run_tomography()


@pytest.mark.parametrize("n_det", [1, 2])
@given(rho=state())
@tomography_settings
def test_json_import_recovers_state(n_det, rho):
    estimate, _, _ = tomography_from_json(json_data(rho, n_det))

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY


@pytest.mark.parametrize("n_det", [1, 2])
@given(rho=state(n_qubits=2), dim=st.floats(0.3, 0.7))
@tomography_settings
def test_json_relative_intensity(n_det, rho, dim):
    """relative_intensity is the drift-correction intensity of a setting (shared by its simultaneous outcomes)."""
    data = json_data(rho, n_det)
    for datum in data["data"]:
        intensity = dim if datum["basis"][0] in MUB_GROUPS[0] else 1.0
        datum["relative_intensity"] = intensity
        datum["counts"][-1] = round(datum["counts"][-1] * intensity)

    estimate, _, fval = tomography_from_json(data)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY
    assert fval < GOOD_FIT


@given(rho=state(n_qubits=2), seed=st.integers(0, 2**31 - 1))
@tomography_settings
def test_json_data_order_does_not_matter(rho, seed):
    data = json_data(rho, 2)
    np.random.default_rng(seed).shuffle(data["data"])

    estimate, _, _ = tomography_from_json(data)

    assert qLib.fidelity(rho, estimate) > EXACT_FIDELITY
