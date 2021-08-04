from TestRun import saveRunsGeneral

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""




"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is used to generate the data for the General analysis. Run this script then run the 
analyze_General.rmd in the results folder"""

# Runs tomography with all setting combinations. Returns the total number of errors occurred.
def runFull(nStates=1,saveStates = True):
    totalErrors = 0
    # I Seperated these two for loops so we don't run invalid settings
    for nBits in [1]:
        for err in [0,3]:
            for acc in [False]:
                for drift in [False,True]:
                    for cross in [False,True]:
                        for bell in [False]:
                            for det in [False,True]:
                                for meth in ["MLE","Linear"]: # Todo add BME
                                    totalErrors += saveRunsGeneral(nBits, nStates,
                                                       resultsFilePath="Results/results_FullData.csv",
                                                       randomStateDist="density",
                                                       errBounds=err, testAccCorr=acc, test2Det=det, testCrossTalk=cross,
                                                       testBell=bell, testDrift=drift, saveData=saveStates, method=meth)
    for nBits in [2]:
        for err in [0,3]:
            for acc in [False,True]:
                for det in [False,True]:
                    for cross in [False,True]:
                        for bell in [False,True]:
                            for drift in [False,True]:
                                for meth in ["MLE","Linear"]: # Todo add BME
                                    totalErrors += saveRunsGeneral(nBits, nStates,
                                                       resultsFilePath="Results/results_FullData.csv",
                                                       randomStateDist="density",
                                                       errBounds=err, testAccCorr=acc, test2Det=det, testCrossTalk=cross,
                                                       testBell=bell, testDrift=drift, saveData=saveStates, method=meth)
    # # Currently accidentals don't work with 3 qubits
    # for nBits in [3]:
    #     for err in [0,3]:
    #         for acc in [False]:
    #             for det in [False,True]:
    #                 for cross in [False,True]:
    #                     for bell in [False,True]:
    #                         for drift in [False,True]:
    #                             for meth in ["MLE","Linear"]: # Todo add BME
    #                                 totalErrors += saveRunsGeneral(nBits, nStates,
    #                                                    resultsFilePath="Results/results_FullData.csv",
    #                                                    randomStateDist="density",
    #                                                    errBounds=err, testAccCorr=acc, test2Det=det, testCrossTalk=cross,
    #                                                    testBell=bell, testDrift=drift, saveData=saveStates, method=meth)

    return totalErrors

if __name__ == '__main__':
    runFull(100)