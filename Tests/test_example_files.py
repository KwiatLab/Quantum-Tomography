import unittest
from pathlib import Path
import numpy as np
import numpy.testing as npt
import QuantumTomography as qLib

"""
Tests that each file in ExampleFiles/ loads, runs, and produces a valid density
matrix.  Where the expected state is known the fidelity is also checked.
"""

EXAMPLES_DIR = Path(__file__).parent.parent / "ExampleFiles"

# Pure state vectors for the idealized examples
_R_STATE = np.array([1, 1j], dtype=complex) / np.sqrt(2)
_PHI_PLUS = np.array([1, 0, 0, 1], dtype=complex) / np.sqrt(2)
_PHI_MINUS = np.array([1, 0, 0, -1], dtype=complex) / np.sqrt(2)


def _assert_valid_density_matrix(rho, tol=1e-6):
    npt.assert_allclose(rho, rho.conj().T, atol=tol, err_msg="rho is not Hermitian")
    npt.assert_allclose(np.trace(rho), 1.0, atol=tol, err_msg="rho trace != 1")
    eigenvalues = np.linalg.eigvalsh(rho)
    assert np.all(eigenvalues >= -tol), f"rho has negative eigenvalue: {eigenvalues.min():.3e}"


def _run(data_file, conf_file=None):
    t = qLib.Tomography()
    if conf_file is not None:
        ext = str(conf_file).split('.')[-1]
        if ext == 'toml':
            t.import_conf(str(conf_file))
        else:
            t.importConf(str(conf_file))

    ext = str(data_file).split('.')[-1]
    if ext == 'json':
        t.import_data(str(data_file))
    else:
        t.importData(str(data_file))

    rho, intensity, fval = t.run_tomography()
    return rho, intensity, fval


class Test_Example_Files(unittest.TestCase):
    def test_1_qubit(self):
        rho, _, _ = _run(EXAMPLES_DIR / "1_qubit_example.json")
        _assert_valid_density_matrix(rho)
        self.assertGreater(qLib.fidelity(_R_STATE, rho), 0.95)

    def test_bell_state(self):
        rho, _, _ = _run(EXAMPLES_DIR / "bell_state_example.json")
        _assert_valid_density_matrix(rho)
        self.assertGreater(qLib.fidelity(_PHI_PLUS, rho), 0.95)

    def test_2_detector(self):
        rho, _, _ = _run(EXAMPLES_DIR / "2n_detector_example.json")
        _assert_valid_density_matrix(rho)
        self.assertGreater(qLib.fidelity(_PHI_PLUS, rho), 0.95)

    def test_crosstalk_and_inefficiency(self):
        rho, _, _ = _run(EXAMPLES_DIR / "crosstalk_inefficiency_example.json")
        _assert_valid_density_matrix(rho)
        self.assertGreater(qLib.fidelity(_PHI_PLUS, rho), 0.90)

    def test_conf_toml_with_data(self):
        # Verifies that a TOML conf file can be loaded alongside a data file
        rho, _, _ = _run(EXAMPLES_DIR / "bell_state_example.json", conf_file=EXAMPLES_DIR / "conf.toml")
        _assert_valid_density_matrix(rho)
        self.assertGreater(qLib.fidelity(_PHI_PLUS, rho), 0.95)

    def test_python_eval(self):
        rho, _, _ = _run(EXAMPLES_DIR/ "pythoneval.txt", conf_file=EXAMPLES_DIR/"conf.txt")
        _assert_valid_density_matrix(rho)
        self.assertLess(abs(qLib.purity(rho) - 1.0), 0.0001) # Arbitray, but to see if this result changes

    def test_python_video_eval(self):
        rho, _, _ = _run(EXAMPLES_DIR/ "pythoneval_video.txt", conf_file=EXAMPLES_DIR/"conf.txt")
        _assert_valid_density_matrix(rho)
        self.assertLess(abs(qLib.purity(rho) - 0.91), 0.01) # Arbitray, but to see if this result changes

    def test_bell_psi_eval(self):
        rho, _, _ = _run(EXAMPLES_DIR/ "bell_psi_data.txt", conf_file=EXAMPLES_DIR/"bell_psi_conf.txt")
        _assert_valid_density_matrix(rho)
        self.assertLess(abs(qLib.purity(rho) - 0.734), 0.003) # Arbitray, but to see if this result changes


if __name__ == "__main__":
    unittest.main()
