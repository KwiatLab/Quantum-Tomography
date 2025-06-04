try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping
import re
from pathlib import Path
import json
import numpy as np

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

DEFAULT_CONF = {
    "n_qubits": 1,
    "n_detectors": 1,
    "crosstalk": np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]),
    "use_derivative": 0,
    "bellstate": 0,
    "do_error_estimation": 0,
    "do_drift_correction": False,
    "do_accidental_correction": False,
    "window": 0,
    "efficiency": np.array([1, 1, 1, 1]),
    "rho_start": None,
    "beta": None,
    "ftol": None,
    "xtol": None,
    "gtol": None,
    "maxfev": None,
    "method": "MLE",
}


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


def export_json(filename, dict):
    filepath = Path(filename)
    with open(filepath) as f:
        json.dump(dict, f)


def import_conf(filename):
    json_dict = import_json(filename)
    # Populate the config with the defaults if they don't exist
    for key, value in DEFAULT_CONF.items():
        if key not in json_dict.keys():
            json_dict[key] = value

    return json_dict


def import_json(filename):
    filepath = Path(filename)
    with open(filepath) as f:
        json_dict = json.load(f)

    for key, value in json_dict.items():
        if isinstance(value, list):
            array = np.asarray(value, dtype=complex)
            json_dict[key] = array

    return json_dict
