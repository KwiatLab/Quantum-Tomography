import numpy as np

import KwiatQuantumLib as qKLib

#create tomography object
t = qKLib.Tomography()


#import data and conf file or just the eval file
t.importConf('ExampleFiles/conf.txt')
[rho, intensity, fval] = t.importData('ExampleFiles/data.txt')
#or just the eval file
[rho, intensity, fval] = t.importPythonEval('ExampleFiles/pythoneval.txt')



#or set the conf settings directly or with the helper
t.conf['DoAccidentalCorrection'] = 1
t.setConfSetting('DoAccidentalCorrection',1)
#and call the object's tomography function
tomo_input=np.array([[1,0,500,1,0],[1,0,0,0,1],[1,0,250,0.7071,0.7071],[1,0,250,0.7071,-0.7071],[1,0,250,0.7071,0.7071j],[1,0,250,0.7071,-0.7071j]])
intensity=np.array([1,1,1,1,1,1])
[rho, intens, fval] = t.state_tomography(tomo_input, intensity)



#The library also include usefull functions you may use. See TomoFunctions.py for more details
expectedState = np.array([[1.0,0.0],[0.0,0.0]])
print("Fidelity: " + str(qKLib.fidelity(rho,expectedState)))
print("State: ")
print(rho)