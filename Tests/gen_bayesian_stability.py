import unittest
from TestRun import runTests
import QuantumTomography as qLib
import matplotlib.pyplot as plt
import numpy as np
import warnings
warnings.filterwarnings("ignore", category=RuntimeWarning)

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""




"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"""This script is UNTRACKED by git, and can be used however. 
By default it returns true. Most likely if you are running tests you'll want to save the
 results and use the SaveRun class, an example of using that class is in the save_RandomStates"""


if __name__ == '__main__':
    numRuns = 25


    maxIter = 0
    numPriorArray = np.zeros(numRuns, dtype=int)
    numPosteArray = np.zeros(numRuns, dtype=int)
    numTimeArray = np.zeros(numRuns, dtype=int)
    numTotalArray = np.zeros(numRuns, dtype=int)

    Tomographys = runTests(1,numRuns, method="BME")
    stabilities = np.zeros((numRuns, len(Tomographys[0][0].stabilityHistory)), dtype=np.double)
    for i in range(numRuns):
        [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = Tomographys[i]
        try:
            numTimeArray[i] = Total_Time
            numPriorArray[i] = Tomo_Object.priorsToConverge
            numPosteArray[i] = Tomo_Object.last_fval - Tomo_Object.priorsToConverge
            numTotalArray[i] = Tomo_Object.last_fval
            stabilities[i] = Tomo_Object.stabilityHistory
            numIterPriors = np.int(np.floor((Tomo_Object.priorsToConverge+1)/1000))
            numIterTotal = np.int(np.floor((Tomo_Object.last_fval + 1) / 1000))
            maxIter = max(maxIter,numIterTotal)

            plt.plot(np.arange(1,1+numIterPriors),
                     np.log10(stabilities[i,1:(numIterPriors+1)]), ',-m', alpha=.25)
            plt.plot(np.arange(numIterPriors,numIterPriors+len(stabilities[i,numIterPriors:numIterTotal])),
                     np.log10(stabilities[i,numIterPriors:numIterTotal]), ',-g', alpha=.25)
            print("Fidelity : " + str(Fidelity_with_Original))

        except:
            pass

    print("")
    print("Avg num from Prior : " + str(np.average(numPriorArray)) + ", sd : " + str(np.std(numPriorArray)))
    print("Avg num from Poste : " + str(np.average(numPosteArray)) + ", sd : " + str(np.std(numPosteArray)))
    print("Avg num from Total : " + str(np.average(numTotalArray)) + ", sd : " + str(np.std(numTotalArray)))
    print("Avg Time : " + str(np.average(numTimeArray)) + ", sd : " + str(np.std(numTimeArray)))

    print("# Tomographies: "+str(numRuns))

    plt.xlabel("Samples (1000)")
    plt.ylabel("Stability log10(y)")
    plt.title("Stability Convergence")
    plt.xlim(1,maxIter)
    plt.show()
