from TestRun import saveRunsGeneral

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""




"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"""This script is used to generate the data for the General analysis. Run this script then run the 
analyze_General.rmd in the results folder"""

if __name__ == '__main__':
    saveRunsGeneral(numQubits=1, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=2, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=3, nStates=100,method='MLE')

    saveRunsGeneral(numQubits=1, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=2, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=1, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=2, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=1, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=2, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=3, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=3, nStates=100,method='MLE')
    saveRunsGeneral(numQubits=3, nStates=100,method='MLE')
    # saveRunsGeneral(numQubits=1, nStates=100,method='BME')
    # saveRunsGeneral(numQubits=2, nStates=100,method='BME')
    # saveRunsGeneral(numQubits=3, nStates=100,method='BME')