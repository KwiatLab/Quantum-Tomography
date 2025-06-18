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
from typing import Union, List, Tuple
from typing_extensions import Self

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


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
    save_state = True

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
    coincidence_window: float = 0.0
    rel_efficiency: List[float] = [0, 0, 0, 0]
    crosstalk: np.ndarray = np.zeros((2**n_qubits, 2**n_qubits))
    measurement_densities: np.ndarray = np.array(
        [density for _, density in POLARIZATION_DENSITIES.items()]
    )
    counts: np.ndarray = np.zeros((n_qubits**n_detectors, n_measurements_per_qubit))
    intensity: np.ndarray = np.ones(len(measurement_densities))
    overall_norms: np.ndarray = np.kron(intensity, rel_efficiency)
    accidentals: np.ndarray = np.zeros_like(counts)

    @model_validator(mode="after")
    def check_input_shape(self) -> Self:
        self.counts = np.atleast_2d(np.array(self.counts, dtype=np.complex128))
        if self.counts.shape != (
            self.n_qubits**self.n_detectors,
            self.n_measurements_per_qubit,
        ):
            raise ValueError(
                f"Counts shape is wrong {self.counts.shape}, expected {self.n_qubits**self.n_detectors, self.n_measurements_per_qubit}"
            )

        return self


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


def import_data(filename: str, config: Union[TomoConfiguration, None] = None):
    filepath = Path(filename)
    with open(filepath) as f:
        json_dict = json.load(f)

    if config:
        json_dict["config"] = config
    data = TomoData.model_validate(json_dict)

    return data


def import_config(filename: str):
    filepath = Path(filename)
    with open(filepath, "rb") as f:
        conf = tomllib.load(f)

    config = TomoConfiguration.model_validate(conf)
    return config
