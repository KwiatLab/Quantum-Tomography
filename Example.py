from __future__ import print_function
from traceback import print_exc
import numpy as np
import matplotlib.pyplot as plt
import QuantumTomography as qtomo

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

np.set_printoptions(precision=3, suppress=True)
# Step 1: Import configuration and data files

tomo_config = qtomo.import_config("ExampleFiles/conf.toml")

tomo_data = qtomo.import_data("ExampleFiles/1_qubit_example.json", tomo_config)

# Step 2: Initialize tomography object

tomo_obj = qtomo.Tomography(tomo_data, tomo_config)

# Step 3. Run Tomography on The data
results = tomo_obj.StateTomography(tomo_data, tomo_config)
print(results)

# import data file
# or call the object's tomography function
tomo_input = np.array(
    [
        [1, 0, 500, 1, 0],
        [1, 0, 0, 0, 1],
        [1, 0, 250, 0.7071, 0.7071],
        [1, 0, 250, 0.7071, -0.7071],
        [1, 0, 250, 0.7071, 0.7071j],
        [1, 0, 250, 0.7071, -0.7071j],
    ]
)
intensity = np.array([1, 1, 1, 1, 1, 1])


# [rho, intens, fval] = t.state_tomography(tomo_input, intensity)
# or import the eval file to import both the config and data
# [rho, intensity, fval] = t.importEval("ExampleFiles/pythoneval.txt")

# Step 4. Optional Methods
# The library also include useful functions you may use included in TomoFunctions.py.
# See https://quantumtomo.web.illinois.edu/Doc/ for a full reference guide of all the functions.
# qKLib.printLastOutput(t)

# expectedState = np.array([[1.0, 0.0], [0.0, 0.0]])
# print("Fidelity: " + str(qKLib.fidelity(rho, expectedState)))
# rho = np.kron(rho, rho)
# qKLib.makeRhoImages(rho, plt, True)
# plt.show()
