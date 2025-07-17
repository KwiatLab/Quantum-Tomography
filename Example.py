from __future__ import print_function
import numpy as np
import matplotlib.pyplot as plt
import QuantumTomography as qKLib

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


# Step 1. Initialize Tomography Object
t = qKLib.Tomography()

# Step 2. Set up Configurations
t.import_conf("ExampleFiles/conf.toml")
t.import_data("ExampleFiles/1_qubit_example.json")

# Step 3. Run Tomography on The data

[rho, intens, fval] = t.run_tomography()

# Step 4. Optional Methods
# The library also include useful functions you may use included in TomoFunctions.py.
# See https://quantumtomo.web.illinois.edu/Doc/ for a full reference guide of all the functions.
qKLib.printLastOutput(t)

expected_state = np.array([[1.0, 0.0], [0.0, 0.0]])
print("Fidelity: " + str(qKLib.fidelity(rho, expected_state)))
qKLib.makeRhoImages(rho, plt, True)
plt.show()


# 2 qubit example

t.import_conf("ExampleFiles/conf.toml")
t.import_data("ExampleFiles/bell_state_example.json")

expected_bell_state = np.array(
    [[0.5, 0, 0, 0.5], [0, 0, 0, 0], [0, 0, 0, 0], [0.5, 0, 0, 0.5]]
)

[rho, intens, fval] = t.run_tomography()

print(rho)

print("Fidelity: " + str(qKLib.fidelity(rho, expected_bell_state)))
qKLib.makeRhoImages(rho, plt, True)
plt.show()


# 2n detector example

t.import_conf("ExampleFiles/conf.toml")
t.import_data("ExampleFiles/2n_detector_example.json")

expected_bell_state = np.array(
    [[0.5, 0, 0, 0.5], [0, 0, 0, 0], [0, 0, 0, 0], [0.5, 0, 0, 0.5]]
)

t.setConfSetting("do_drift_correction", True)
[rho, intens, fval] = t.run_tomography()


print("Fidelity: " + str(qKLib.fidelity(rho, expected_bell_state)))
qKLib.makeRhoImages(rho, plt, True)
plt.show()

###################################################################
### LEGACY CODE BELOW, WILL NOT BE SUPPORTED IN FUTURE VERSIONS ###
###################################################################


# If you still have a need to use the old input file formats, you can do so with the same import functions as before.
# Note that this will be removed in future versions.

t = qKLib.Tomography()

# Step 2. Set up Configurations
# import conf file
t.importConf('ExampleFiles/conf.txt')
# or set the conf settings directly 
t.conf["DoAccidentalCorrection"] = 0

# Step 3. Run Tomography on The data

# import data file
# importing the data file will automatically run the tomography
[rho, intensity, fval] = t.importData('ExampleFiles/data.txt')

qKLib.makeRhoImages(rho, plt, True)
plt.show()

[rho,intensity,fval] = t.importEval('ExampleFiles/pythoneval.txt')
qKLib.makeRhoImages(rho, plt, True)
plt.show()

