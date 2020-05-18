import numpy as np

import KwiatQuantumLib as qKLib

# Step 1. Initialize Tomography Object
t = qKLib.Tomography()

#Step 2. Set up Configurations
# import conf file
t.importConf('ExampleFiles/conf.txt')
#or set the conf settings directly or with the helper
t.setConfSetting('DoAccidentalCorrection',1)

# Step 3. Run Tomography on The data
# import data file
[rho, intensity, fval] = t.importData('ExampleFiles/data.txt')
# or call the object's tomography function
tomo_input=np.array([[1,0,500,1,0],[1,0,0,0,1],[1,0,250,0.7071,0.7071],[1,0,250,0.7071,-0.7071],[1,0,250,0.7071,0.7071j],[1,0,250,0.7071,-0.7071j]])
intensity=np.array([1,1,1,1,1,1])
[rho, intens, fval] = t.state_tomography(tomo_input, intensity)
# or import the eval file to import both the config and data
[rho, intensity, fval] = t.importEval('ExampleFiles/pythoneval.txt')


# Step 4. Optional Methods
# The library also include useful functions you may use. See TomoFunctions.py for more details
expectedState = np.array([[1.0,0.0],[0.0,0.0]])
print("Fidelity: " + str(qKLib.fidelity(rho,expectedState)))
print("State: ")
print(rho)