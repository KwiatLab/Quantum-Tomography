"""Regression tests for bugs fixed in 1.2.2."""

import copy
import json
from pathlib import Path

import matplotlib
import numpy as np
import numpy.testing as npt
import pytest

import QuantumTomography as qLib

matplotlib.use("Agg")
EXAMPLES_DIR = Path(__file__).parent.parent / "ExampleFiles"
BELL = json.loads((EXAMPLES_DIR / "bell_state_example.json").read_text())
TWO_DET = json.loads((EXAMPLES_DIR / "2n_detector_example.json").read_text())


def noiseless_counts(rho, measurements):
    """Counts for a 1 detector per qubit measurement of rho."""
    n = measurements.shape[1] // 2
    kets = [np.array(1)] * len(measurements)
    for k in range(n):
        kets = [np.kron(ket, row[2 * k : 2 * k + 2]) for ket, row in zip(kets, measurements)]
    return np.array([1e4 * np.real(ket.conj() @ rho @ ket) for ket in kets])


def run_json(data):
    t = qLib.Tomography()
    t._import_data(copy.deepcopy(data))
    return t, t.run_tomography()[0]


@pytest.mark.parametrize("n_qubits", [1, 2])
def test_linear_returns_state_not_transpose(n_qubits):
    t = qLib.Tomography(n_qubits)
    rho = qLib.random_density_state(n_qubits)
    measurements = t.getStandardBasis(n_qubits, 1)
    rho_linear = t.StateTomography(measurements, noiseless_counts(rho, measurements), method="LINEAR")[0]
    npt.assert_allclose(rho_linear, rho, atol=1e-6)


@pytest.mark.parametrize("method", ["LINEAR", "MLE"])
def test_global_phase_of_measurement_kets_is_irrelevant(method):
    t = qLib.Tomography(2)
    rho = qLib.random_density_state(2)
    measurements = t.getStandardBasis(2, 1)
    counts = noiseless_counts(rho, measurements)
    phases = np.exp(2j * np.pi * np.random.rand(len(measurements), 2))
    rho_phased = t.StateTomography(measurements * np.repeat(phases, 2, axis=1), counts, method=method)[0]
    assert qLib.fidelity(rho_phased, rho) > 0.999


def test_2det_relative_intensity_is_applied():
    scaled = copy.deepcopy(TWO_DET)
    for datum in scaled["data"]:
        brighter = set(datum["basis"]) <= {"H", "V"}
        datum["relative_intensity"] = 3.0 if brighter else 1.0
        datum["counts"] = [c * datum["relative_intensity"] for c in datum["counts"]]
    npt.assert_allclose(run_json(scaled)[1], run_json(TWO_DET)[1], atol=1e-4)


def test_2det_data_order_is_irrelevant():
    shuffled = copy.deepcopy(TWO_DET)
    np.random.shuffle(shuffled["data"])
    npt.assert_allclose(run_json(shuffled)[1], run_json(TWO_DET)[1], atol=1e-4)


def test_2det_orthogonal_partner_found_with_rounding_error():
    data = copy.deepcopy(TWO_DET)
    data["measurement_states"]["D"] = [np.cos(0.3), np.sin(0.3)]
    data["measurement_states"]["A"] = [np.sin(0.3), -np.cos(0.3)]
    run_json(data)


def test_3_qubit_2det_counts_are_recognised():
    t = qLib.Tomography(3)
    measurements = t.getStandardBasis(3, 2)
    t.buildTomoInput(measurements, np.ones((len(measurements), 8)), -1, 0, -1, -1, 0, 0)
    assert t.conf["NDetectors"] == 2


def test_1det_json_with_relative_intensity_and_no_window():
    data = copy.deepcopy(BELL)
    del data["coincidence_window"]
    for datum in data["data"]:
        datum["relative_intensity"] = 1.0
    npt.assert_allclose(run_json(data)[1], run_json(BELL)[1], atol=1e-6)


def test_accidental_correction_3_qubits():
    t = qLib.Tomography(3)
    t.conf["DoAccidentalCorrection"] = 1
    measurements = t.getStandardBasis(3, 1)
    n = len(measurements)
    t.intensities = np.ones(n)
    t.last_input = t.buildTomoInput(
        measurements, np.ones(n), -1, 0, np.ones(n), np.full((n, 3), 1e3), np.array([1.0]), 0
    )
    accidentals = t.filter_data(t.last_input)[3]
    npt.assert_allclose(accidentals, 1e3**3 * 1e-9**2)


def test_wrong_length_time_raises():
    t = qLib.Tomography(1)
    with pytest.raises(ValueError, match="Invalid time array"):
        t.StateTomography(t.getStandardBasis(1, 1), np.ones(6), time=np.ones(3))


def test_bell_settings_error_bars():
    t, _ = run_json(BELL)
    kwargs = dict(partsize_init=3, partsize=3, t=1)
    names, means, errors = t.getBellSettings(bounds=3, **kwargs).T
    expected = np.array(
        [qLib.getBellSettings_helper(state[0], **kwargs)[:, 1] for state in t.mont_carlo_states], dtype=float
    )
    assert list(names) == ["s", "a", "a'", "b", "b'"]
    npt.assert_allclose(means.astype(float), expected.mean(axis=0))
    npt.assert_allclose(errors.astype(float), expected.std(axis=0, ddof=1))


def test_negativity_of_separable_states_is_zero():
    assert qLib.negativity(np.eye(4) / 4) == 0
    assert qLib.negativity(qLib.toDensity(np.array([1, 0, 0, 0], dtype=complex))) == 0


def test_remove_global_phase_normalises():
    npt.assert_allclose(qLib.removeGlobalPhase(np.array([2j, 0])), [1, 0], atol=1e-12)


def test_conf_accepts_n_for_no():
    t = qLib.Tomography()
    t.conf["Bellstate"] = "n"
    assert t.conf["Bellstate"] == 0


def test_rho_images_3_qubits():
    import matplotlib.pyplot as plt

    qLib.makeRhoImages(np.eye(8) / 8, plt)
    plt.close("all")


def test_linear_subtracts_accidentals():
    t = qLib.Tomography(2)
    t.conf["DoAccidentalCorrection"] = 1
    rho = qLib.random_density_state(2)
    measurements = t.getStandardBasis(2, 1)
    n = len(measurements)
    accidentals = 1e5**2 * 100e-9
    rho_linear = t.StateTomography(
        measurements,
        noiseless_counts(rho, measurements) + accidentals,
        time=np.ones(n),
        singles=np.full((n, 2), 1e5),
        window=np.array([100.0]),
        method="LINEAR",
    )[0]
    npt.assert_allclose(rho_linear, rho, atol=1e-6)


def test_2det_json_without_detectors_used_raises():
    data = copy.deepcopy(TWO_DET)
    for datum in data["data"]:
        del datum["detectors_used"]
    with pytest.raises(ValueError, match="detectors_used"):
        run_json(data)


def test_hedge_does_not_depend_on_intensity():
    t_vals = qLib.density2t(qLib.random_density_state(1))
    measurements = np.array([np.eye(2)], dtype=complex)

    def fitness(scale):
        counts = np.array([scale**2])  # equals the prediction, so only the hedge term is left
        return qLib.maxlike_fitness_hedged(scale * t_vals, counts, np.zeros(1), measurements, 0.5, np.ones(1))

    npt.assert_allclose(fitness(1), fitness(100))


def test_purity_of_ket():
    npt.assert_allclose(qLib.purity(np.array([1, 1j]) / np.sqrt(2)), 1)


@pytest.mark.parametrize("n", [0, 1, 2])
def test_partial_transpose_of_each_subsystem(n):
    parts = [qLib.random_density_state(1) for _ in range(3)]
    expected = [part.T if k == n else part for k, part in enumerate(parts)]
    kron3 = lambda a, b, c: np.kron(np.kron(a, b), c)
    npt.assert_allclose(qLib.partial_transpose(kron3(*parts), n), kron3(*expected))
