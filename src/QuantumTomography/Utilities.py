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


class TomoConfiguration(BaseModel):
    n_qubits: int = 1
    n_detectors: int = 1
    n_measurements_per_qubit: int = 6
    use_derivative: bool = False
    get_bell_settings: bool = False
    coincidence_window: float = 0.0
    rel_efficiency: List[float] = [0, 0, 0, 0]
    do_error_estimation: bool = False
    do_drift_correction: bool = False
    do_accidental_correction: bool = False
    beta: Union[float, None] = None
    ftol: Union[float, None] = None
    xtol: Union[float, None] = None
    gtol: Union[float, None] = None
    maxfev: Union[int, None] = None
    method: TomographyType = TomographyType.MLE


class TomoData(BaseModel):
    class Config:
        arbitrary_types_allowed = True

    config: TomoConfiguration = TomoConfiguration()
    measurement_densities: np.ndarray = np.array(
        [density for _, density in POLARIZATION_DENSITIES.items()]
    )
    counts: Union[List, np.ndarray] = np.zeros(
        (config.n_qubits**config.n_detectors, config.n_measurements_per_qubit)
    )
    crosstalk: np.ndarray = np.zeros((2**config.n_qubits, 2**config.n_qubits))

    @model_validator(mode="after")
    def check_input_shape(self) -> Self:
        config = self.config
        self.counts = np.atleast_2d(np.array(self.counts, dtype=np.complex128))
        if self.counts.shape != (
            config.n_qubits**config.n_detectors,
            config.n_measurements_per_qubit,
        ):
            raise ValueError(
                f"Counts shape is wrong {self.counts.shape}, expected {config.n_qubits**config.n_detectors, config.n_measurements_per_qubit}"
            )

        return self


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


class ConfDict(MutableMapping):
    """A dictionary where the casing of the keys don't matter
    and values can be automatically handled."""

    def __init__(self, *args, **kwargs):
        self.store = dict()
        self.update(dict(*args, **kwargs))  # use the free update to set keys

    def __getitem__(self, key):
        return self.store[self._keytransform(key)]

    def __setitem__(self, key, value):
        self.store[self._keytransform(key)] = self._valuetransform(value)

    def __delitem__(self, key):
        del self.store[self._keytransform(key)]

    def __iter__(self):
        return iter(self.store)

    def __len__(self):
        return len(self.store)

    def _keytransform(self, key):
        return key.lower()

    def _valuetransform(self, value):
        if isinstance(value, str):
            if (
                value.lower() == "yes"
                or value.lower() == "true"
                or value.lower() == "t"
                or value.lower() == "y"
            ):
                value = 1
            elif (
                value.lower() == "no"
                or value.lower() == "false"
                or value.lower() == "f"
                or value.lower() == "f"
            ):
                value = 0
            elif (
                value.upper() == "LINEAR"
                or value.upper() == "MLE"
                or value.upper() == "HMLE"
                or value.upper() == "BME"
            ):
                return value.upper()
            else:
                raise ValueError('Invalid Conf Setting of "' + value + '"')

        return value


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
