try:
    from collections.abc import MutableMapping
except ImportError:
    from collections import MutableMapping
import re

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



def getValidFileName(fileName):
    newFileName = re.sub(r"^\\.+", "", fileName)
    newFileName = re.sub(r"[\\\\/:*?\"<>|]", "", newFileName)
    if newFileName == "":
        raise ValueError("File Name : '" + fileName + "' results in an empty fileName!");
    return newFileName