from TestRun import runTests
import QuantumTomography as qLib
import os
import numpy as np

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is UNTRACKED by git, and can be used however. 
By default it returns true. Most likely if you are running tests you'll want to save the
 results and use the SaveRun class, an example of using that class is in the save_RandomStates"""

if __name__ == '__main__':

    states = [qLib.random_pure_state(1),qLib.random_bell_state(1),qLib.random_density_state(1),
              qLib.random_pure_state(2),qLib.random_bell_state(2),qLib.random_density_state(2)]

    resultsFilePath = "Results/results_Monte_Carlo.csv"

    # Open a file with access mode 'a'
    file_object = open(resultsFilePath, 'w')
    printHeader=True
    # Do tomographies of single qubit
    for state in states:
        [[Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time]] = runTests(int(np.floor(np.log2(state.shape[0])+.01)),1,
                                                                                        randomStateDist=state)
        # todo: try error corrections with 2det and acc
        for i in [10,9,8,7,6,5,4,3,2,1,0]:
            props = Tomo_Object.getProperties(i)
            dataRow = dict()
            dataRow['num_errors'] = str(i)
            dataRow['intensity_mean'] = props[0,1]
            dataRow['intensity_std'] = props[0,2]
            dataRow['concurrence_mean'] = props[2,1]
            dataRow['concurrence_std'] = props[2,2]
            dataRow['entropy_mean'] = props[4,1]
            dataRow['entropy_std'] = props[4,2]
            dataRow['negativity_mean'] = props[6,1]
            dataRow['negativity_std'] = props[6,2]
            dataRow['purity_mean'] = props[7,1]
            dataRow['purity_std'] = props[7,2]
            if printHeader:
                TORREPLACE = ""
                for key in dataRow.keys():
                    TORREPLACE += key + ","
                printHeader = False
            else:
                TORREPLACE=","

            TORREPLACE= TORREPLACE[:-1] +"\n"
            for key in dataRow.keys():
                TORREPLACE += str(dataRow[key])+","
            TORREPLACE = TORREPLACE[:-1]

            # Append at the end of file
            file_object.write(TORREPLACE)

    # Close the file
    file_object.close()