try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping
import re
from pathlib import Path
from enum import Enum, auto
import json
import tomllib
import numpy as np
from pydantic import BaseModel, model_validator, field_validator
from typing import Union, List, Tuple, get_type_hints, Type, Dict, Any, Optional
from typing_extensions import Self
import functools
from math import comb

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

POLARIZATION_STATE_NAMES = ["H", "V", "D", "A", "R", "L"]


class Polarizations(Enum):
    H = auto()
    V = auto()
    D = auto()
    A = auto()
    R = auto()
    L = auto()


POLARIZATION_STATES = {
    Polarizations.H: np.array([1, 0], dtype=np.complex128),
    Polarizations.V: np.array([0, 1], dtype=np.complex128),
    Polarizations.D: 1 / np.sqrt(2) * np.array([1, 1], dtype=np.complex128),
    Polarizations.A: 1 / np.sqrt(2) * np.array([1, -1], dtype=np.complex128),
    Polarizations.R: 1 / np.sqrt(2) * np.array([1, 1.0j], dtype=np.complex128),
    Polarizations.L: 1 / np.sqrt(2) * np.array([1, -1.0j], dtype=np.complex128),
}

POLARIZATION_DENSITIES = {}
for name, state in POLARIZATION_STATES.items():
    POLARIZATION_DENSITIES[name] = np.outer(state, state.conj().T)

print(POLARIZATION_DENSITIES)


class TomographyType(Enum):
    LINEAR = 1
    MLE = 2
    HMLE = 3
    BME = 4


TOMOGRAPHY_TYPES = {
    "LINEAR": TomographyType.LINEAR,
    "MLE": TomographyType.MLE,
    "HMLE": TomographyType.HMLE,
    "BME": TomographyType.BME,
}


class TomoConfiguration(BaseModel):
    """Model describing how the tomography should be performed. This is used along with TomoData. This is loaded from TOML.

    Fields:

    use_derivative (bool):

    get_bell_settings (bool):

    do_error_estimation (bool):

    do_drift_correction (bool):

    do_accidental_correction (bool):

    beta (float, optional):

    ftol (float, optional):
    xtol (float, optional):
    gtol (float, optional):
    maxfev (int, optional):
    method (TomographyType or str):
    starting_matrix (np.ndarray, optional):
    save_state (bool)
    minimizer_kwargs (dict):
    """

    class Config:
        arbitrary_types_allowed = True

    use_derivative: bool = False
    get_bell_settings: bool = False
    do_error_estimation: bool = False
    do_drift_correction: bool = False
    do_accidental_correction: bool = False
    beta: Optional[float] = None
    ftol: Optional[float] = None
    xtol: Optional[float] = None
    gtol: Optional[float] = None
    maxfev: Optional[int] = None
    method: Union[TomographyType, str] = TomographyType.MLE
    starting_matrix: Union[np.ndarray, None] = None
    save_state: bool = True
    _minimizer_kwargs: Dict = {
        "ftol": ftol,
        "xtol": xtol,
        "gtol": gtol,
        "maxfev": maxfev,
    }

    @model_validator(mode="after")
    def check_method(self) -> Self:
        if isinstance(self.method, str):
            if self.method not in TOMOGRAPHY_TYPES.keys():
                raise ValueError(
                    f"Tomography Type {self.method} not supported! Must be one of {[key for key in TOMOGRAPHY_TYPES.keys()]}"
                )
            else:
                self.method = TOMOGRAPHY_TYPES[self.method]
        return self

    @model_validator(mode="after")
    def get_not_none_minimizer_kwargs(self) -> Self:
        self._minimizer_kwargs = {
            k: v for k, v in self._minimizer_kwargs.items() if v is not None
        }
        return self


def get_highest_fold_coincidence_count_index(n_fold):
    idx = 0
    for i in range(n_fold):
        idx += comb(n_fold, i)
    return idx - 1


class Measurement(BaseModel):
    """Model describing individual tomographic measurements.

    Fields:

    basis (List[str]): The set of projector used for a measurement.
        For multiple qubits, this will be a list longer than 1.
        Ex: [H,V] measures H on qubit 1 and V on qubit 2 simultaneously.

    integration_time: The integration time of the detectors used for this measurement. Used for accidental correction.

    counts: The singles and coincidence counts for the detectors.
        For a multi-qubit measurement,the first len(basis) entries are the singles and the rest are coincidences.
        Ex 1: basis = [H,V], counts = [{H counts}, {V counts}, {H-V coincidences}]
        Ex 2: basis = [H,V,D], counts = [{H counts}, {V counts}, {H-V coincidences}, {H-D coincidences}, {H-V-D coincidences}]

    accidentals: A square 2d array of size (n_detectors_per_qubit**n_qubits, n_detectors_per_qubit**n_qubits) describing the accidental counts between each pair of detectors.
        Because of this, the matrix is symmetric and the diagonal ignored. Note that this is only used for coincidence measurements

    detectors_used: A list of indices describing what detectors are used for this specific measurement. This is used for accidental correction, in correlation
        to the coincidence_window field in TomoData. It is also used for calculation of relative detector pair inefficiency in correlation with relative_efficiency.
        Note: this is not used when there is only 1 detector per measurement.

    relative_intensity: Relative intensity of this measurement. Used to correct for intensity drift. Note that this is not used when n_detectors_per_qubit > 1. This scales the counts relative to the other measurements.
    """

    class Config:
        arbitrary_types_allowed = True

    basis: List[str]
    integration_time: float

    counts: np.ndarray

    overall_norm: Optional[float] = 1.0

    accidentals: int = 0

    relative_intensity: Optional[float] = 1.0

    detectors_used: List[int] = [0]

    _crosstalk_corrected_density: Optional[np.ndarray] = None

    _coincidence_counts: Optional[np.ndarray] = None

    @field_validator("counts", mode="before")
    @classmethod
    def cast_to_ndarray(cls, counts: Any) -> np.ndarray:
        return np.array(counts)

    @model_validator(mode="after")
    def get_coincidences_only(self) -> Self:
        self._coincidence_counts = np.array(
            self.counts[get_highest_fold_coincidence_count_index(len(self.basis))]
        )
        return self


class TomoData(BaseModel):
    """Model describing a collection of tomographic measurements. This is loaded from JSON.

    Fields:

    config: A TomoConfiguration object, describing how the tomography should be run on this data.

    n_qubits: Number of qubits in this dataset.

    n_detectors: Number of detectors used for measurements in this dataset.

    relative_efficiency: A square 2d array of size (n_detectors**n_qubits, n_detectors**n_qubits) describing the relative efficiencies between each pair of detectors.
        Because of this, the matrix will be symmetric and the diagonal ignored.
        Note that this is only used when n_detectors > 1 for the purpose of coincidence inefficiency correction.

    crosstalk:  A square 2d array of size (n_measurements_per_qubit * n_qubits, n_measurements_per_qubit * n_qubits) describing the crosstalk between each measurement.
        This matrix may NOT be column stochastic to allow for absorbance correction.
        Note the convention of this matrix is dependent on the order of the data given.
        Ex. if data = ["HH", "HV", "VH", "VV"],
            then crosstalk[1,2] is the probability that measurement "HV" will be measured as "VH"

    coincidence_window: A square 2d array of size (n_detectors**n_qubits, n_detectors**qubits) describing the coincidence window between each pair of detectors.

    measurement_densities: A dictionary mapping names (strings) of measurements to projectors (ndarray).
        This allows users to name and define their own measurement projectors used in their experiment.

    simultaneous_measurement_indices: A list of lists of indices describing which measurements are done simultaneously. Note this is only used for n_detectors > 1.
        This is used to normalize counts between orthogonal and simultaneous measurements to correct for intensity drift.
        Ex:
            measurement_densities = {"H", "V", "D", "A", "R","L"},
            simultaneous_measurement_indices = [[0,1], [2,3], [4,5]]

    data: A list of Measurement objects.
    """

    class Config:
        arbitrary_types_allowed = True

    config: TomoConfiguration = TomoConfiguration()
    n_qubits: int = 1
    n_detectors_per_qubit: int = 1

    n_measurements_per_qubit: int = 6

    relative_efficiency: np.ndarray = np.array([1])

    crosstalk: np.ndarray = np.eye(
        n_measurements_per_qubit * n_qubits, n_measurements_per_qubit * n_qubits
    )

    coincidence_window: np.ndarray = np.zeros(
        (n_qubits * n_detectors_per_qubit, n_qubits * n_detectors_per_qubit)
    )

    measurement_densities: Dict[str, np.ndarray] = {
        POLARIZATION_STATE_NAMES[i]: POLARIZATION_DENSITIES[name]
        for i, name in enumerate(POLARIZATION_STATES)
    }

    simultaneous_measurement_indices: Optional[List[List[int]]] = None

    _all_measurement_densities: Optional[np.ndarray] = None

    data: List[Measurement] = [
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[0]],
            integration_time=1.0,
            counts=np.array([50]),
            accidentals=0,
            relative_intensity=1.0,
            detectors_used=[0],
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[1]],
            integration_time=1.0,
            counts=np.array([50]),
            accidentals=0,
            relative_intensity=1.0,
            detectors_used=[0],
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[2]],
            integration_time=1.0,
            counts=np.array([50]),
            accidentals=0,
            relative_intensity=1.0,
            detectors_used=[0],
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[3]],
            integration_time=1.0,
            counts=np.array([50]),
            accidentals=0,
            relative_intensity=1.0,
            detectors_used=[0],
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[4]],
            integration_time=1.0,
            counts=np.array([100]),
            accidentals=0,
            relative_intensity=1.0,
            detectors_used=[0],
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[5]],
            integration_time=1.0,
            counts=np.array([0]),
            accidentals=0,
            relative_intensity=1.0,
            detectors_used=[0],
        ),
    ]

    @model_validator(mode="after")
    def calculate_overall_norm(self) -> Self:
        if self.config.do_drift_correction:
            for datum in self.data:
                datum.overall_norm = datum.relative_intensity
                if self.n_detectors_per_qubit > 1:
                    det_1, det_2 = datum.detectors_used
                    datum.overall_norm *= self.relative_efficiency[det_1, det_2]

        return self

    @field_validator("coincidence_window", mode="after")
    @classmethod
    def make_coincidence_window_symmetric(cls, coincidence_window) -> np.ndarray:
        coincidence_window = coincidence_window + coincidence_window.T
        return coincidence_window

    @field_validator("relative_efficiency", mode="after")
    @classmethod
    def make_rel_efficiency_symmetric(cls, relative_efficiency) -> np.ndarray:
        relative_efficiency = relative_efficiency + relative_efficiency.T
        return relative_efficiency

    @model_validator(mode="after")
    def calculate_crosstalk_corrected_densities(self) -> Self:
        """Get corrected measurement projectors that account for crosstalk.

        For more info see Altpeter, J. et al., "Photonic State Tomography", p.33-34.
        https://research.physics.illinois.edu/QI/Photonics/tomography-files/amo_tomo_chapter.pdf
        """

        meas_basis_len = len(self.data[0].basis)
        density_shape = (2**meas_basis_len, 2**meas_basis_len)

        # Precompute the measurement densities using the basis definitions
        all_densities = get_all_densities_from_data(self)

        # Do crosstalk correction
        for i in range(meas_basis_len):
            corrected_density = np.zeros(density_shape, dtype=np.complex128)
            for j in range(meas_basis_len):
                # New densities are linear combinations of all measurement projectors multiplied by a scalar
                corrected_density += self.crosstalk[i, j] * all_densities[j]

            self.data[i]._crosstalk_corrected_density = corrected_density

        return self

    @field_validator("data", mode="before")
    @classmethod
    def cast_data_to_measurement(cls, value: List[dict]) -> List[Measurement]:
        new_list = []
        for datum in value:
            if isinstance(datum, dict):
                new_list.append(Measurement(**datum))
            else:
                new_list.append(datum)
        return new_list

    @model_validator(mode="after")
    def cast_lists_to_array(self) -> Self:
        for key, density in self.measurement_densities.items():
            self.measurement_densities[key] = np.array(density)

        self.measurement_densities = self.measurement_densities
        self.relative_efficiency = np.array(self.relative_efficiency)
        return self

    @model_validator(mode="after")
    def calculate_accidentals(self) -> Self:
        """Calculate accidentals to use for correction.
        For more info see Altpeter, J. et al., "Photonic State Tomography", p.33.
        https://research.physics.illinois.edu/QI/Photonics/tomography-files/amo_tomo_chapter.pdf
        """

        if self.config.do_accidental_correction == 1:
            for datum in self.data:
                # multiply all of the singles together
                accidental_count_factor = np.prod(datum.counts[0 : len(datum.basis)])
                # Multiply by the coincidence window and divide by integration time
                datum.accidentals = (
                    accidental_count_factor
                    * self.coincidence_window[
                        datum.detectors_used[0], datum.detectors_used[1]
                    ]
                    / datum.integration_time
                )

        return self


def get_all_densities_from_data(tomo_data: TomoData) -> np.ndarray:
    all_densities = []
    for datum in tomo_data.data:
        # Get all of the densities used in this Measurement
        densities = [tomo_data.measurement_densities[name] for name in datum.basis]
        # Kronecker product all of them together
        all_densities.append(functools.reduce(lambda x, y: np.kron(x, y), densities))

    return np.array(all_densities)


def cast_to_numpy(json_dict, key):
    json_dict[key] = np.array(json_dict[key], dtype=np.complex128)


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


def getValidFileName(fileName):
    newFileName = re.sub(r"^\\.+", "", fileName)
    newFileName = re.sub(r"[\\\\/:*?\"<>|]", "", newFileName)
    if newFileName == "":
        raise ValueError("File Name : '" + fileName + "' results in an empty fileName!")
    return newFileName


def export_data(filename, tomo_data: TomoData):
    filepath = Path(filename)
    with open(filepath) as f:
        json.dump(dict, f)


def get_fields_annotations(m: Type[BaseModel]) -> Dict[str, Any]:
    default_annotations = get_type_hints(m)
    return {
        field_name: default_annotations[field_name] for field_name in m.model_fields
    }


def import_data(filename: str, config: TomoConfiguration):
    filepath = Path(filename)
    with open(filepath) as f:
        json_dict = json.load(f)

    for density_key, density_value in json_dict["measurement_densities"].items():
        cast_to_numpy(json_dict["measurement_densities"], density_key)

    data_fields = get_fields_annotations(TomoData)
    for key, value in data_fields.items():
        print(value)
        if value is np.ndarray and key in json_dict:
            print("casting")
            cast_to_numpy(json_dict, key)

    data = TomoData(**json_dict, config=config)
    print(data)
    return data


def import_config(filename: str):
    filepath = Path(filename)
    with open(filepath, "rb") as f:
        conf = tomllib.load(f)

    config = TomoConfiguration.model_validate(conf)
    return config
