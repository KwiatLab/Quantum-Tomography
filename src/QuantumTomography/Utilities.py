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
from pydantic import BaseModel, model_validator
from typing import Union, List, Tuple, get_type_hints, Type, Dict, Any
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
    class Config:
        arbitrary_types_allowed = True

    use_derivative: bool = False
    get_bell_settings: bool = False
    do_error_estimation: bool = False
    do_drift_correction: bool = False
    do_accidental_correction: bool = False
    beta: Union[float, None] = None
    ftol: Union[float, None] = None
    xtol: Union[float, None] = None
    gtol: Union[float, None] = None
    maxfev: Union[int, None] = None
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


class TomoData(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    config: TomoConfiguration = TomoConfiguration()
    n_qubits: int = 1
    n_detectors: int = 1
    n_measurements_per_qubit: int = 6
    rel_efficiency: np.ndarray = np.array([1])
    crosstalk: np.ndarray = np.zeros((2**n_qubits, 2**n_qubits))

    coincidence_window: np.ndarray = np.zeros((n_detectors, n_detectors))
    measurement_densities: Dict[str, np.ndarray] = {
        POLARIZATION_STATE_NAMES[i]: POLARIZATION_DENSITIES[name]
        for i, name in enumerate(POLARIZATION_STATES)
    }
    intensity: np.ndarray = np.ones(len(measurement_densities))
    overall_norms: np.ndarray = np.kron(rel_efficiency, intensity)

    data: List[Dict] = [
        {
            "basis": [POLARIZATION_STATE_NAMES[0]],
            "integration_time": 1,
            "counts": [50],
            "accidentals": np.zeros(
                len(np.choose(np.arange(0, n_detectors), np.arange(0, n_detectors)))
            ),
        },
        {
            "basis": [POLARIZATION_STATE_NAMES[1]],
            "integration_time": 1,
            "counts": [50],
            "accidentals": np.zeros(
                len(np.choose(np.arange(0, n_detectors), np.arange(0, n_detectors)))
            ),
        },
        {
            "basis": [POLARIZATION_STATE_NAMES[2]],
            "integration_time": 1,
            "counts": [50],
            "accidentals": np.zeros(
                len(np.choose(np.arange(0, n_detectors), np.arange(0, n_detectors)))
            ),
        },
        {
            "basis": [POLARIZATION_STATE_NAMES[3]],
            "integration_time": 1,
            "counts": [50],
            "accidentals": np.zeros(
                len(np.choose(np.arange(0, n_detectors), np.arange(0, n_detectors)))
            ),
        },
        {
            "basis": [POLARIZATION_STATE_NAMES[4]],
            "integration_time": 1,
            "counts": [100],
            "accidentals": np.zeros(
                len(np.choose(np.arange(0, n_detectors), np.arange(0, n_detectors)))
            ),
        },
        {
            "basis": [POLARIZATION_STATE_NAMES[5]],
            "integration_time": 1,
            "counts": [0],
            "accidentals": np.zeros(
                len(np.choose(np.arange(0, n_detectors), np.arange(0, n_detectors)))
            ),
        },
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

    @model_validator(mode="after")
    def cast_lists_to_array(self) -> Self:
        for key, density in self.measurement_densities.items():
            self.measurement_densities[key] = np.array(density)

        self.measurement_densities = np.array(self.measurement_densities)
        self.rel_efficiency = np.array(self.rel_efficiency)
        for datum in self.data:
            datum["data"] = np.array(datum["data"])
            datum["accidentals"] = np.array(datum["accidentals"])
        return self

    @model_validator(mode="after")
    def check_accidental_shape(self) -> Self:
        if self.config.do_accidental_correction == 1:
            expected_accidental_shape = len(
                np.choose(np.arange(0, n_detectors), np.arange(0, n_detectors))
            )
            for datum in self.data:
                if datum["accidentals"].shape != expected_accidental_shape:
                    raise ValueError(
                        f"Accidentals aren't the correct shape. Expected {expected_accidental_shape} but got {datum['accidentals'].shape} for basis {datum['basis']}."
                    )

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


def import_data(filename: str, config: Union[TomoConfiguration, None] = None):
    filepath = Path(filename)
    with open(filepath) as f:
        json_dict = json.load(f)

    data_fields = get_fields_annotations(TomoData)
    for key, value in data_fields.items():
        print(value)
        if value is np.ndarray and key in json_dict:
            print("casting")
            cast_to_numpy(json_dict, key)

    data = TomoData(**json_dict)
    print(data)
    return data


def import_config(filename: str):
    filepath = Path(filename)
    with open(filepath, "rb") as f:
        conf = tomllib.load(f)

    config = TomoConfiguration.model_validate(conf)
    return config
