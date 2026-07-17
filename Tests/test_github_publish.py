from __future__ import print_function
import warnings

import pytest
from gen_Full import runFull

warnings.filterwarnings("ignore")


"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"""This test script is the one used by github when a new version is created. This will automatically
run and the results can be see in the actions tab"""

"Attention! These tests run on the version that your environment uses. see readme for details"


#    1 Qubit
def test_FULL():
    # TODO: Write better test
    pytest.skip("Disabled: runFull is slow and not run as part of the regular suite.")
    numErrors = runFull(nStates=1, saveStates=False)
    assert numErrors == 0
