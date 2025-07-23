try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping
import re
import numpy as np
import functools
from math import comb
from pathlib import Path
import ast

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


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
        if (isinstance(value, str)):
            if (value.lower() == "yes" or
                value.lower() == "true" or
                value.lower() == "t" or
                value.lower() == "y"):
                value = 1
            elif (value.lower() == "no" or
                  value.lower() == "false" or
                  value.lower() == "f"or
                  value.lower() == "f"):
                value = 0
            elif(value.upper() == "LINEAR" or
                value.upper() == "MLE" or
                value.upper() == "HMLE" or
                value.upper() == "BME"):
                return value.upper()
            else:
                raise ValueError('Invalid Conf Setting of "' + value+'"')

        return value
    
def parse_np_array(string):
    array = string.strip().replace("(","").replace(")","").strip("np.array")
    array = array.strip(";")
    # Stolen directly from stack overflow https://stackoverflow.com/a/43879517
    return np.array(ast.literal_eval(re.sub(r'\]\s*\[', r'],[', re.sub(r'(\d+)\s+(\d+)', r'\1,\2', array.replace('\n','')))))

def get_raw_measurement_bases_from_data(tomo_data) -> np.ndarray:
    all_densities = []
    for datum in tomo_data["data"]:
        # Get all of the densities used in this Measurement
        densities = [
            np.array(tomo_data["measurement_states"][name], dtype=np.complex128)
            for name in datum["basis"]
        ]
        all_densities.append(np.array(densities).flatten())

    return np.array(all_densities)

def get_all_product_states_from_data(tomo_data) -> np.ndarray:
    all_densities = []
    for datum in tomo_data["data"]:
        # Get all of the densities used in this Measurement
        densities = [
            np.array(tomo_data["measurement_states"][name], dtype=np.complex128)
            for name in datum["basis"]
        ]
        # Kronecker product all of them together
        all_densities.append(functools.reduce(lambda x, y: np.kron(x, y), densities))

    return np.array(all_densities)


def get_highest_fold_coincidence_count_index(n_fold):
    idx = 0
    for i in range(n_fold):
        idx += comb(n_fold, i)
    return idx - 1


def cast_to_numpy(json_dict, key):
    json_dict[key] = np.array(json_dict[key], dtype=np.complex128)


def getValidFileName(fileName):
    newFileName = re.sub(r"^\\.+", "", fileName)
    newFileName = re.sub(r"[\\\\/:*?\"<>|]", "", newFileName)
    if newFileName == "":
        raise ValueError("File Name : '" + fileName + "' results in an empty fileName!");
    return newFileName
