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

expected_state = np.array([1,1.j])
expected_state = np.outer(expected_state, expected_state.conj())
print("Fidelity: " + str(qKLib.fidelity(rho, expected_state)))
qKLib.makeRhoImages(rho, plt, True)
plt.show()


# Data import/export example
t = qKLib.Tomography()

# Set up Configurations
t.import_conf("ExampleFiles/conf.toml")
t.import_data("ExampleFiles/1_qubit_example.json")

# You can export the configurations from the new format (toml/json) to the old format (.txt)
t.exportToData("ExampleFiles/data_exported.txt")
t.exportToConf("ExampleFiles/conf_exported.txt")
t.exportToEval("ExampleFiles/eval_exported.txt")

# You can also export the configuration to the new format for later use
t.export_to_conf_toml("ExampleFiles/conf_toml_exported.toml")

t2 = qKLib.Tomography()
t2.importConf("conf_exported.txt")
t2.importData("data_exported.txt")
[rho, intens, fval] = t2.run_tomography()

qKLib.printLastOutput(t2)

expected_state = np.array([1,1.j])
expected_state = np.outer(expected_state, expected_state.conj())
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

# example with crosstalk and detector pair inefficiencies

t.import_conf("ExampleFiles/conf.toml")
t.import_data("ExampleFiles/crosstalk_inefficiency_example.json")


t.conf["DoDriftCorrection"] = 1
expected_bell_state = np.array(
    [[0.5, 0, 0, 0.5], [0, 0, 0, 0], [0, 0, 0, 0], [0.5, 0, 0, 0.5]]
)

[rho, intens, fval] = t.run_tomography()


print(f"Optimal CHSH bell measurement settings: {t.getBellSettings()}")

print("Fidelity: " + str(qKLib.fidelity(rho, expected_bell_state)))
qKLib.makeRhoImages(rho, plt, True)
plt.show()


###################################################################
### LEGACY CODE BELOW, WILL NOT BE SUPPORTED IN FUTURE VERSIONS ###
###################################################################


# If you still have a need to use the old input file formats, you can do so with the same import functions as before.
# Note that this will be removed in future versions.

# t = qKLib.Tomography()
#
## Step 2. Set up Configurations
## import conf file
# t.importConf('ExampleFiles/conf.txt')
## or set the conf settings directly
# t.conf["DoAccidentalCorrection"] = 0
#
## Step 3. Run Tomography on The data
#
## import data file
## importing the data file will automatically run the tomography
# [rho, intensity, fval] = t.importData('ExampleFiles/data.txt')
#
# qKLib.makeRhoImages(rho, plt, True)
# plt.show()
#
# [rho,intensity,fval] = t.importEval('ExampleFiles/pythoneval.txt')
# qKLib.makeRhoImages(rho, plt, True)
# plt.show()
