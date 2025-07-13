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
    """Model describing how the tomography should be performed. This is used along with TomoData.

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
    minimizer_kwargs: Dict = {
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

    accidentals: A square 2d array of size (n_detectors**n_qubits, n_detectors**n_qubits) describing the accidental counts between each pair of detectors.
        Because of this, the matrix is symmetric and the diagonal ignored. Note that this is only used for coincidence measurements

    relative_intensity: Relative intensity of this measurement. Used to correct for intensity drift. Note that this is not used when n_detectors > 1. This scales the counts relative to the other measurements.
    """

    class Config:
        arbitrary_types_allowed = True

    basis: List[str]
    integration_time: float

    counts: Union[np.ndarray, List[int]]

    accidentals: Optional[np.ndarray] = None

    relative_intensity: Optional[float] = 1.0

    @field_validator("counts", mode="before")
    @classmethod
    def cast_to_ndarray(cls, counts: Any) -> np.ndarray:
        return np.array(counts)


class TomoData(BaseModel):
    """Model describing a collection of tomographic measurements

    Fields:

    config: A TomoConfiguration object, describing how the tomography should be run on this data.

    n_qubits: Number of qubits in this dataset.

    n_detectors: Number of detectors used for measurements in this dataset.

    rel_efficiency: A square 2d array of size (n_detectors**n_qubits, n_detectors**n_qubits) describing the relative efficiencies between each pair of detectors.
        Because of this, the matrix will be symmetric and the diagonal ignored.
        Note that this is only used when n_detectors > 1 for the purpose of coincidence inefficiency correction.

    crosstalk:  A square 2d array of size (n_detectors**n_qubits, n_detectors**n_qubits) describing the crosstalk between each detector.
        This matrix may NOT be column stochastic to allow for absorbance correction.

    coincidence_window: A square 2d array of size (n_detectors**n_qubits, n_detectors**qubits) describing the coincidence window between each pair of detectors.

    measurement_densities: A dictionary mapping names (strings) of measurements to projectors (ndarray).
        This allows users to name and define their own measurement projectors used in their experiment.

    orthogonal_measurement_indices: A list of tuples describing which measurements are done simultaneously. Note this is only used for n_detectors > 1.
        This is used to normalize counts between orthogonal and simultaneous measurements to correct for intensity drift.
        Ex:
            measurement_densities = {"H", "V", "D", "A", "R","L"},
            orthogonal_measurement_indices = [(0,1), (2,3), (4,5)]

    data: A list of Measurement objects.
    """

    class Config:
        arbitrary_types_allowed = True

    config: TomoConfiguration = TomoConfiguration()
    n_qubits: int = 1
    n_detectors: int = 1

    n_measurements_per_qubit: int = 6

    rel_efficiency: np.ndarray = np.array([1])

    crosstalk: np.ndarray = np.eye(n_qubits * n_detectors, n_qubits * n_detectors)

    coincidence_window: np.ndarray = np.zeros(
        (n_qubits * n_detectors, n_qubits * n_detectors)
    )

    measurement_densities: Dict[str, np.ndarray] = {
        POLARIZATION_STATE_NAMES[i]: POLARIZATION_DENSITIES[name]
        for i, name in enumerate(POLARIZATION_STATES)
    }

    orthogonal_measurement_indices: Optional[List[Tuple[int]]] = None

    data: List[Measurement] = [
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[0]],
            integration_time=1.0,
            counts=[50],
            accidentals=np.zeros((n_detectors**n_qubits, n_detectors**n_qubits)),
            relative_intensity=1.0,
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[1]],
            integration_time=1.0,
            counts=[50],
            accidentals=np.zeros((n_detectors**n_qubits, n_detectors**n_qubits)),
            relative_intensity=1.0,
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[2]],
            integration_time=1.0,
            counts=[50],
            accidentals=np.zeros((n_detectors**n_qubits, n_detectors**n_qubits)),
            relative_intensity=1.0,
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[3]],
            integration_time=1.0,
            counts=[50],
            accidentals=np.zeros((n_detectors**n_qubits, n_detectors**n_qubits)),
            relative_intensity=1.0,
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[4]],
            integration_time=1.0,
            counts=[100],
            accidentals=np.zeros((n_detectors**n_qubits, n_detectors**n_qubits)),
            relative_intensity=1.0,
        ),
        Measurement(
            basis=[POLARIZATION_STATE_NAMES[5]],
            integration_time=1.0,
            counts=[0],
            accidentals=np.zeros((n_detectors**n_qubits, n_detectors**n_qubits)),
            relative_intensity=1.0,
        ),
    ]
    # @model_validator(mode="after")
    # def check_input_shape(self) -> Self:
    # self.counts = np.atleast_2d(np.array(self.counts, dtype=np.complex128))
    # if self.counts.shape != (
    #    self.n_qubits**self.n_detectors,
    #    self.n_measurements_per_qubit,
    # ):
    #    raise ValueError(
    #        f"Counts shape is wrong {self.counts.shape}, expected {self.n_qubits**self.n_detectors, self.n_measurements_per_qubit}"
    #    )

    #    return self

    # @model_validator(mode="after")
    # def make_window_symmetric(self)->Self:

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
            print("casting11")
            self.measurement_densities[key] = np.array(density)

        self.measurement_densities = self.measurement_densities
        self.rel_efficiency = np.array(self.rel_efficiency)
        return self

    @model_validator(mode="after")
    def check_accidental_shape(self) -> Self:
        if self.config.do_accidental_correction == 1:
            expected_accidental_shape = len(
                np.choose(
                    np.arange(0, self.n_detectors), np.arange(0, self.n_detectors)
                )
            )
            for datum in self.data:
                if datum.accidentals.shape != expected_accidental_shape:
                    raise ValueError(
                        f"Accidentals aren't the correct shape. Expected {expected_accidental_shape} but got {datum['accidentals'].shape} for basis {datum['basis']}."
                    )
        return self

    @model_validator(mode="after")
    def accidental_correction(self) -> Self:
        if self.config.do_accidental_correction == 1:
            scalerIndex = np.concatenate((np.ones(self.n_qubits - 2), [2, 2]))
            additiveIndex = np.array([0, 1])
            for j in range(2, self.n_qubits):
                additiveIndex = np.concatenate(([2 * j], additiveIndex))
            if len(coinc.shape) == 1:
                acc = acc[:, np.newaxis]
            for j in range(n_coinc):
                index = bin(j).split("b")[1]
                index = "0" * (nbits - len(index)) + index
                index = [int(char) for char in index]
                index = index * scalerIndex + additiveIndex
                index = np.array(index, dtype=int)
                acc[:, j] = np.prod(np.real(sings[:, tuple(index)]), axis=1) * (
                    window[j] * 1e-9 / np.real(times)
                ) ** (nbits - 1)
            if acc.shape != coinc.shape:
                acc = acc[:, 0]
        return self


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
