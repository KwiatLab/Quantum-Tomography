# Quantum-Tomography
![Python package](https://github.com/KwiatQIM/Quantum-Tomography/workflows/Python%20package/badge.svg?branch=master)

A python library to help perform tomography on a quantum state.

Links
 - [Documentation](http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/)
 - [Video Tutorial](https://www.youtube.com/watch?v=I-214P0LOfQ&list=PLJLHMKtk5Pqy9w9aCuyowUF1p7pl2JCI9&index=3)
 - [Pypi Package](https://pypi.org/project/Quantum-Tomography/)


## Usage
### Terminal
For those who do not want to write python code to perform tomography 
on their data, you can use the following command in the package directory:
```
Quantum-Tomography -i C:\Full\Path\To\pythoneval.txt
```
This will read the data in the txt file provided and print the output ot the console. Examples and syntax for a conf file
is provided at the bottom of this readme. If you would like to save your data you can provide the save location like in the following example:
```
Quantum-Tomography -i C:\Full\Path\To\pythoneval.txt -s C:\Full\Path\To\output\Folder 
```
There are several other arguments that are optional like save option. Here is the full list of arguments:

- -i or --eval
    - Values : string
    - Desc : The full path to the file that contains the data and configuration for the tomography.
- -s or --save
    - Values : string
    - Desc : The full path to the folder where you want the output to be saved. If not included it will not save your data.
- -p or --pic
    - Desc : Including this will show images of real and imaginary values of the density matrix. If save is also included pictures will only be saved and not shown.
    - Default : False
    
### Python
For those running tomography on multiple quantum states it may be easier to use the python 
package directly for a more hands free process.
## Usage
##### Step 1. Initialize Tomography Object
```
import QuantumTomography as qKLib

t = qKLib.Tomography()
```
##### Step 2. Set up Configurations
This can be done in multiple ways. The first and the easiest is to pass the full path
to importConf Function. Examples and syntax for a conf file is provided at the bottom of
this readme.
```
t.importConf('Path/To/conf.txt')
```
Specific Config settings can also be set directly using the following function:
```
t.setConfSetting('DoAccidentalCorrection',1)
```
Please refrain from setting the conf settings directly, this does not handle invalid arguments and
can lead to errors. A list values for config is provided at the bottom of this readme and
also in the TomoClass.py file.
##### Step 3. Run Tomography on The data
This can also be done in multiple ways. The first is using the importData Function. Examples
and syntax for a data file is provided at the bottom of this readme and also in the TomoClass.py file..
```
[rho, intens, fval] = t.importData('Path/To/data.txt')
```
Data settings can also be passed into the main tomography function
```
tomo_input = np.array([[1,0,500,1,0],[1,0,0,0,1],[1,0,250,0.7071,0.7071],[1,0,250,0.7071,-0.7071],[1,0,250,0.7071,0.7071j],[1,0,250,0.7071,-0.7071j]])
intensity = np.array([1,1,1,1,1,1])
[rho, intens, fval] = t.state_tomography(tomo_input, intensity)
```
Steps 2 and 3 can be done in one single step by passing in a eval file.
```
[rho, intensity, fval] = t.importEval('Path/To/pythoneval.txt')
```
For running multiple states with the same settings it is recommended to run the tomographying using
the python eval method since the the configurations is being unnecessarily being reset every time.
Examples and syntax for a eval file is provided at the bottom of this readme.
##### Step 4. Optional Methods
We provide many functions to help describe the state. Properties of the state can be calculated with methods found in
the TomoFunctions files. Some examples of these are proveded
```
fid = qKLib.fidelity(rho,expectedState)
entangle = qKLib.entanglement(rho)
entropy = qKLib.entropy(rho)
```
---
## Syntax
### Conf File
This file states the configurations of the tomography. The syntax of the txt file is python. You write the
conf settings just like you would set a python dictionary. These are the following values you can set in a conf file.

- 'NQubits'
    - Values : >= 1
    - Desc : The number of qubits the quantum state has. It will take exponentially more time for more qubits.
    - Default : 2
- 'NDetectors'
    - Values : 1 or 2
    - Desc : The number of detectors per qubit used during the physical tomography of the quantum state.
    - Default : 1
- 'ctalk'
    - Values : matrix that is (2^NQubits) by (2^NQubits)
    - Desc : Cross talk Matrix of the setup.
    - Default : identity matrix with appropriate size
- 'Bellstate'
    - Values : 'no' or 'yes'
    - Desc : Give the optimal measurement settings for a CHSH bell inequality for the estimated density matrix.
    These settings are found through a numerical search over all possible measurement settings.
    - Default : 'no'
- 'DoDriftCorrection'
    - Values : 'no' or 'yes'
    - Desc : Whether of not you want to perform drift correction on the state
    - Default : 'no'
- 'DoAccidentalCorrection'
    - Values : 'no' or 'yes'
    - Desc : Whether of not you want to perform accidental corrections on the state.
    - Default : 'no'
- 'DoErrorEstimation'
    - Values : >=0
    - Desc : Number of states used to calculate the errors on the properties of the state
    - Default : 0
- 'Window'
    - Values : 0 or array like, dimension = 1
    - Desc : Coincidence window durations (in nanoseconds) to calculate the accidental rates. The
    four windows should be entered in the order of the detector pairs 1-2, 1-4, 3-2, 3-4, where A-B
    corresponds to a coincidence measurement between detector A and detector B.
    - Default : '0'
- 'Efficiency'
    - Values : 0 or array like, dimension = 1
    - Desc :  vector that lists the relative coincidence efficiencies of detector pairs when using 2 detectors per
    qubit. The order is detector 1-2, 1-4, 3-2, 3-4.
    - Default : 0
##### Example:
```
conf['NQubits'] = 2
conf['NDetectors'] = 1
conf['Crosstalk'] = [[0.9842,0.0049,0.0049,0],[0.0079,0.9871,0,0.0050],[0.0079,0,0.9871,0.0050],[0.001,0.0079,0.0079,0.9901]]
conf['UseDerivative'] = 0
conf['Bellstate'] = 1
conf['DoErrorEstimation'] = 3
conf['DoDriftCorrection'] = 'no'
conf['Window'] = 0
conf['Efficiency'] = [0.9998,1.0146,0.9195,0.9265]
```
### Data File
This file states the data of the measurements. Both tomo_input the intensity must be specified. The syntax of the txt file is python. You write the
data settings just like you would set a python matrix. This is the following layout of the tomo_input matrix
- tomo_input
    - Values : numpy array, dimension = 2
    - Desc : Relative pump power (arb. units) during measurement; used for drift correction.
    #### For n detectors:
    - tomo_input[:, 0]: times
    - tomo_input[:, 1 : n_qubit + 1)]: singles
    - tomo_input[:, n_qubit + 1]: coincidences
    - tomo_input[:, n_qubit + 2 : 3 * n_qubit + 2)]: measurements

    #### For 2n detectors:
    - tomo_input[:, 0]: times
    - tomo_input[:, 1 : 2 * n_qubit + 1]: singles
    - tomo_input[:, 2 * n_qubit+1 : 2 ** n_qubit + 2 * n_qubit + 1]: coincidences
    - tomo_input[:, 2 ** n_qubit + 2 * n_qubit + 1 : 2 ** n_qubit + 4 * n_qubit + 1 ]: measurements
- intensity
    - Values : numpy array
    - Desc : Relative pump power (arb. units) during measurement; used for drift correction.

##### Example:
This example is for 2 qubits using 1 detector.
```
tomo_input = np.array([[1,0,0,3708,1,0,1,0],[1,0,0,77,1,0,0,1],[1,0,0,1791,1,0,0.7071,0.7071],[1,0,0,2048,1,0,0.7071,0.7071j],[1,0,0,51,0,1,1,0],[1,0,0,3642,0,1,0,1],[1,0,0,2096,0,1,0.7071,0.7071],[1,0,0,1926,0,1,0.7071,0.7071j],[1,0,0,1766,0.7071,0.7071,1,0],[1,0,0,1914,0.7071,0.7071,0,1],[1,0,0,1713,0.7071,0.7071,0.7071,0.7071],[1,0,0,3729,0.7071,0.7071,0.7071,0.7071j],[1,0,0,2017,0.7071,0.7071j,1,0],[1,0,0,1709,0.7071,0.7071j,0,1],[1,0,0,3686,0.7071,0.7071j,0.7071,0.7071],[1,0,0,2404,0.7071,0.7071j,0.7071,0.7071j]])
intensity = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
```
### Eval File
This text file contains all the information for a tomography. It is essentially a conf file and a data file combined into one file.
##### Example:
```
tomo_input = np.array([[1,0,0,3708,1,0,1,0],[1,0,0,77,1,0,0,1],[1,0,0,1791,1,0,0.7071,0.7071],[1,0,0,2048,1,0,0.7071,0.7071j],[1,0,0,51,0,1,1,0],[1,0,0,3642,0,1,0,1],[1,0,0,2096,0,1,0.7071,0.7071],[1,0,0,1926,0,1,0.7071,0.7071j],[1,0,0,1766,0.7071,0.7071,1,0],[1,0,0,1914,0.7071,0.7071,0,1],[1,0,0,1713,0.7071,0.7071,0.7071,0.7071],[1,0,0,3729,0.7071,0.7071,0.7071,0.7071j],[1,0,0,2017,0.7071,0.7071j,1,0],[1,0,0,1709,0.7071,0.7071j,0,1],[1,0,0,3686,0.7071,0.7071j,0.7071,0.7071],[1,0,0,2404,0.7071,0.7071j,0.7071,0.7071j]])
intensity = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
conf['NQubits'] = 2
conf['NDetectors'] = 1
conf['Crosstalk'] = [[0.9842,0.0049,0.0049,0],[0.0079,0.9871,0,0.0050],[0.0079,0,0.9871,0.0050],[0.001,0.0079,0.0079,0.9901]]
conf['UseDerivative'] = 0
conf['Bellstate'] = 0
conf['DoErrorEstimation'] = 1
conf['DoDriftCorrection'] = 'no'
conf['Window'] = 0
conf['Efficiency'] = [0.9998,1.0146,0.9195,0.9265]
```
