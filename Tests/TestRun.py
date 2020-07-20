from __future__ import print_function
import numpy as np
import sys
import matplotlib.pyplot as plt
sys.path.append("../")
import QuantumTomography as qLib
from TestResultsClass import TestResult
import os
import traceback

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""


def runTest(args):
    [numQubits, errBounds, testAccCorr, test2Det, testCrossTalk, testBell, testDrift, nStates] = args
    success = False
    try:
        # set up the test
        tomo = qLib.Tomography()
        allTests = np.zeros(nStates, dtype="O")
        AtotalCounts = np.zeros(nStates)
        myFidels = np.zeros(nStates)

        # set up settings for tomo class
        intensity = np.ones(6 ** numQubits)
        tomo.conf['NQubits'] = numQubits
        tomo.conf['Properties'] = ['concurrence', 'tangle', 'entanglement', 'entropy', 'linear_entropy', 'negativity']
        if (testAccCorr):
            tomo.conf['DoAccidentalCorrection'] = 1
        else:
            tomo.conf['DoAccidentalCorrection'] = 0
        if (test2Det):
            tomo.conf['NDetectors'] = 2
        else:
            tomo.conf['NDetectors'] = 1
        if (not testCrossTalk):
            tomo.conf['Crosstalk'] = np.identity(2 ** numQubits)
        tomo.conf['UseDerivative'] = 0
        tomo.conf['Bellstate'] = testBell
        tomo.conf['DoErrorEstimation'] = errBounds
        if (testDrift):
            tomo.conf['DoDriftCorrection'] = 1
        else:
            tomo.conf['DoDriftCorrection'] = 0
        tomo.conf['Window'] = 1
        tomo.conf['Efficiency'] = np.ones(2 ** numQubits)

        tomo_input = tomo.getTomoInputTemplate()

        # set up measurements
        # measurements is an array of all the measurements for both classes mStates is only to help calculate Measurements
        if (test2Det):
            measurements = np.zeros((len(tomo_input), 2 ** numQubits, 2 ** (numQubits)), dtype=complex)
        else:
            measurements = np.zeros((len(tomo_input), 2 ** numQubits), dtype=complex)
        if (tomo.getNumDetPerQubit() == 1):
            # input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
            mStates = tomo_input[:, np.arange(numQubits + 2, 3 * numQubits + 2)]
        else:
            # input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements
            mStates = tomo_input[:, np.arange(2 ** numQubits + 2 * numQubits + 1,
                                              2 ** numQubits + 4 * numQubits + 1)]
        mStates = np.reshape(mStates, (mStates.shape[0], int(mStates.shape[1] / 2), 2))
        if (test2Det):
            wavePlateArraysBasis = np.zeros((3, 2, 2), dtype=complex)
            wavePlateArraysBasis[0] = np.identity(2, dtype=complex)
            wavePlateArraysBasis[1] = np.array([[.7071, .7071], [.7071, -.7071]], dtype=complex)
            wavePlateArraysBasis[2] = np.array([[.7071, -.7071j], [.7071, .7071j]], dtype=complex)
            wavePlateArray = np.array([1], dtype=complex)
        else:
            wavePlateArraysBasis = np.zeros((6, 2, 2), dtype=complex)
            wavePlateArraysBasis[0] = np.identity(2, dtype=complex)
            wavePlateArraysBasis[1] = np.array([[0, 1], [1, 0]], dtype=complex)
            wavePlateArraysBasis[2] = np.array([[.7071, .7071], [.7071, -.7071]], dtype=complex)
            wavePlateArraysBasis[3] = np.array([[.7071, -.7071], [.7071, .7071]], dtype=complex)
            wavePlateArraysBasis[4] = np.array([[.7071, -.7071j], [.7071, .7071j]], dtype=complex)
            wavePlateArraysBasis[5] = np.array([[.7071, .7071j], [.7071, -.7071j]], dtype=complex)
            wavePlateArray = np.array([1], dtype=complex)
        for i in range(numQubits):
            wavePlateArray = np.kron(wavePlateArray, wavePlateArraysBasis)

        for i in range(len(mStates)):
            if (test2Det):
                for x in range(0, 2 ** (numQubits)):
                    if (x == 1):
                        try:
                            mStates[i][1] = getOppositeState(mStates[i][1])
                        except:
                            mStates[i][0] = getOppositeState(mStates[i][0])
                    elif (x == 2):
                        mStates[i][1] = getOppositeState(mStates[i][1])
                        mStates[i][0] = getOppositeState(mStates[i][0])
                    elif (x == 3):
                        mStates[i][1] = getOppositeState(mStates[i][1])
                    temp = mStates[i][0]
                    for j in range(1, len(mStates[i])):
                        temp = np.kron(temp, mStates[i][j])
                    measurements[i][x] = temp
                try:
                    mStates[i][1] = getOppositeState(mStates[i][1])
                except:
                    mStates[i][0] = getOppositeState(mStates[i][0])
                mStates[i][0] = getOppositeState(mStates[i][0])
            else:
                temp = mStates[i][0]
                for j in range(1, len(mStates[i])):
                    temp = np.kron(temp, mStates[i][j])
                measurements[i] = temp
    except:
        FAIL = '\033[91m'
        print(f"{FAIL}Failed to set up Test: "+uniqueID(args))
        raise SetUpError()

    numErrors = 0
    # create states and do tomo
    for x in range(nStates):
        # crostalk
        cTalkMat = 0
        if (testCrossTalk):
            cTalkMat = np.random.rand(2 ** numQubits, 2 ** numQubits)
            for i in range(2 ** numQubits):
                cTalkMat[:, i] = cTalkMat[:, i] / (sum(cTalkMat[:, i] + 0 * np.random.random()))
            tomo.conf['Crosstalk'] = cTalkMat

        # counts
        numCounts = int(np.random.random() * 10 ** np.random.randint(4, 5))

        # create random state
        state = np.random.beta(.5, .5, 2 ** numQubits) + 1j * np.random.beta(.5, .5, 2 ** numQubits)
        state = state / np.sqrt(np.dot(state, state.conj()))

        startingRho = qLib.toDensity(state)
        # Testing setting
        for i in range(len(tomo_input)):
            # state goes through wave plates
            newState = wavePlateArray[i] @ state
            # state goes through beam splitter and we measure the H counts
            if (testCrossTalk):
                newState = cTalkMat @ newState
            hBasis = np.zeros(2 ** numQubits, dtype=complex)
            hBasis[0] = 1
            if (test2Det):
                prob = np.zeros(2 ** numQubits, complex)
                for j in range(0, 2 * numQubits):
                    hBasis = np.zeros(2 ** numQubits, dtype=complex)
                    hBasis[j] = 1
                    prob[j] = np.dot(hBasis, newState)
                    prob[j] = min(prob[j] * prob[j].conj(), .99999999)
                prob = np.array(prob, dtype=float)
                tomo_input[i,
                2 * numQubits + 1: 2 ** numQubits + 2 * numQubits + 1] = np.random.multinomial(
                    numCounts, prob)
            else:
                prob = np.dot(hBasis, newState)
                prob = prob * prob.conj()
                tomo_input[i, numQubits + 1] = np.random.binomial(numCounts, min(prob, .99999999))
            # input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences

        if (testAccCorr):
            # acc[:, j] = np.prod(np.real(sings[:, idx]), axis=1) * (window[j] * 1e-9 / np.real(t)) ** (nbits - 1)
            if (test2Det):
                sings = tomo_input[:, np.arange(1, 2 * numQubits + 1)]
                coinc = tomo_input[:,
                        np.arange(2 * numQubits + 1, 2 ** numQubits + 2 * numQubits + 1)]
                window = np.ones(4) * 100000
                n_coinc = 2 ** numQubits
            else:
                sings = tomo_input[:, np.arange(1, numQubits + 1)]
                coinc = tomo_input[:, numQubits + 1]
                window = [100000]
                n_coinc = 1
            sings = 2000 * np.random.random(sings.shape)
            t = np.random.random(tomo_input.shape[0]) + 1
            acc = np.zeros_like(coinc)
            if (len(acc.shape) == 1):
                acc = acc[:, np.newaxis]

            scalerIndex = np.concatenate((np.ones(numQubits - 2), [2, 2]))
            additiveIndex = np.array([0, 1])
            for j in range(2, numQubits):
                additiveIndex = np.concatenate(([2 * j], additiveIndex))
            for j in range(n_coinc):
                index = bin(j).split("b")[1]
                index = "0" * (numQubits - len(index)) + index
                index = [int(char) for char in index]
                index = index * scalerIndex + additiveIndex
                index = np.array(index, dtype=int)
                acc[:, j] = np.prod(np.real(sings[:, tuple(index)]), axis=1) * (window[j] * 1e-9 / np.real(t)) ** (
                            numQubits - 1)
            if (acc.shape != coinc.shape):
                acc = acc[:, 0]
            if (test2Det):
                tomo_input[:, np.arange(1, 2 * numQubits + 1)] = sings
                tomo_input[:,
                np.arange(2 * numQubits + 1, 2 ** numQubits + 2 * numQubits + 1)] = coinc + acc
            else:
                tomo_input[:, np.arange(1, numQubits + 1)] = sings
                tomo_input[:, numQubits + 1] = coinc + acc
            tomo.conf['Window'] = window
            tomo_input[:, 0] = t

        if (testDrift):
            intensity = np.random.random(intensity.shape) ** 2 * 10
            if (test2Det):
                # tomo_input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
                for k in range(tomo_input.shape[0]):
                    tomo_input[
                        k, np.arange(2 * numQubits + 1, 2 ** numQubits + 2 * numQubits + 1)] = int(
                        (intensity[k]) * tomo_input[
                            k, np.arange(2 * numQubits + 1, 2 ** numQubits + 2 * numQubits + 1)])
            else:
                # tomo_input[:, n_qubit+1]: coincidences
                for k in range(tomo_input.shape[0]):
                    tomo_input[k, numQubits + 1] = int((intensity[k]) * tomo_input[k, numQubits + 1])
        errorMessage = ""
        # Do tomography with settings
        try:
            myDensitie, inten, myfVal = tomo.state_tomography(tomo_input, intensity)

            myFidel = qLib.fidelity(startingRho, myDensitie)
            if (testBell):
                tomo.getBellSettings(myDensitie)
                tomo.getProperties(myDensitie)
        except:
            myDensitie = [[0]]
            inten = 0
            myfVal = -1
            myFidel = -1
            errorMessage = "<pre>\n" + traceback.format_exc() + "\n</pre>"

        allTests[x] = TestResult(startingRho, myDensitie, myfVal, myFidel, inten, numCounts, cTalkMat, errorMessage)

        AtotalCounts[x] = numCounts
        myFidels[x] = myFidel
        if (myFidel < .8):
            numErrors += 1

    if(numErrors>0):
        FAIL = '\033[91m'
        print(f"{FAIL}At least 1 tomography failed with settings:"+uniqueID(args))
        print(f"{FAIL}View at: ErrorLogs/" + uniqueID(args) + '.html')

        # Print settings used
        FORREPLACE = '<table>'
        FORREPLACE += '<tr><th>Settings used:</th><th></th></tr>'
        FORREPLACE += '<tr><td><ul>'
        settingsArray = [testAccCorr, test2Det, testCrossTalk, testDrift]
        settingsArrayName = ['testAccCorr', 'test2Det', 'testCrossTalk', 'testDrift']
        for x in range(len(settingsArray)):
            if settingsArray[x]:
                FORREPLACE += '<li>'
                FORREPLACE += settingsArrayName[x]
                FORREPLACE += '</li>'
        FORREPLACE += '</ul></td></tr></table>'

        # Print Data
        FORREPLACE += '<table class="data">'
        # HEADERS
        FORREPLACE += '<tr>'
        FORREPLACE += '<th colspan="2">Actual Densities</th><th>Counts</th>'  # Actual
        FORREPLACE += '<th colspan="3">Caclulated Densities</th>'  # Calculated
        if (testCrossTalk):  # Ctalk
            FORREPLACE += '<th colspan="2">CrossTalk</th>'
        FORREPLACE += '</tr>'
        for i in allTests:
            FORREPLACE += str(i)
        FORREPLACE += '</table>'

        with open('ErrorLogs/Template.html', 'r') as f:
            fff = '\n'.join(f.readlines())
            f.close()
        fff = fff.replace('TOREPLACE', str(FORREPLACE))
        with open('ErrorLogs/' + uniqueID(args) + '.html', 'w') as ff:
            ff.write(fff)
            ff.close()

        raise TomographyError()


def uniqueID(args):
    [numQubits, errBounds, testAccCorr, test2Det, testCrossTalk, testBell, testDrift, nStates] = args

    s = "N"+str(numQubits) +"-"
    s +="e"+str(errBounds) +"-"
    if (testAccCorr):
        s+="a1-"
    else:
        s += "a0-"
    if (test2Det):
        s+="d1-"
    else:
        s += "d0-"
    if (testCrossTalk):
        s+="c1-"
    else:
        s += "c0-"
    if (testBell):
        s+="b1-"
    else:
        s += "b0-"
    if (testDrift):
        s+="dr1"
    else:
        s += "dr0"
    return s


class SetUpError(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return "Error Setting up a test"

class TomographyError(Exception):
    def __init__(self):
        pass
    def __str__(self):
        return "Error running tomography on a data set"

def getOppositeState(psi):
    # Horizontal
    if (all(psi == np.array([1, 0], dtype=complex))):
        return np.array([0, 1], dtype=complex)
    if (all(psi == np.array([0, 1], dtype=complex))):
        return np.array([1, 0], dtype=complex)
    # Diagional
    if (all(psi == np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2))], dtype=complex))):
        return np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2))], dtype=complex)
    if (all(psi == np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2))], dtype=complex))):
        return np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2))], dtype=complex)
    # Circle
    if (all(psi == np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j], dtype=complex))):
        return np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j], dtype=complex)
    if (all(psi == np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j], dtype=complex))):
        return np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j], dtype=complex)
    else:
        raise Exception('State Not Found getOppositeState')
