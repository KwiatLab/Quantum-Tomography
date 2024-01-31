from __future__ import print_function
import numpy as np
import quantum_tomography as qLib
import traceback
import warnings
import time
import os

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"""This script runs a multiple tomographies and does not saves the results"""
"""This script is used by test_github_publish and _update."""

"Attention! These tests run on the version that your environment uses. see readme for details"


"""
runTests
This function runs a number of tomography's and returns the objects as well as other things such as the original
state and total time it took to run the tomography.

Parameters
----------
randomStateDist : ["density","pure","bellstate"] or a ndarray
    What distribution you want to use to create the random state. You can also give it a state to use
    
Returns
-------
BigListOfTomographies: list
    Each row has the following form [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time]

"""
def runTests(numQubits, nStates,
             randomStateDist="density",errBounds=0,testAccCorr=False,test2Det=False,testCrossTalk=False,
             testBell=False,testDrift=False,method='MLE'):
    BigListOfTomographies = list()

    # If invalid conf settings then don't run
    if numQubits != 2:
        if testAccCorr:
            print('Run skipped. Invalid Conf settings code: n1-a1')
            return BigListOfTomographies
        if testBell:
            print('Run skipped. Invalid Conf settings code: n'+str(numQubits)+'-b1')
            return BigListOfTomographies

    numErrors = 0
    for TomographyNumber in range(nStates):

        ################################
        # SETUP THE CLASS AND SETTINGS #
        ################################

        # Create Class
        Tomo_Object = qLib.Tomography()
        # set up settings for tomography class
        Tomo_Object.conf['NQubits'] = numQubits
        if (testAccCorr):
            Tomo_Object.conf['DoAccidentalCorrection'] = 'yEs'
        else:
            Tomo_Object.conf['DoAccidentalCorrection'] = "nO"
        if (test2Det):
            Tomo_Object.conf['NDetectors'] = 2
            Tomo_Object.conf['Efficiency'] = np.ones(2 ** numQubits)
        else:
            Tomo_Object.conf['NDetectors'] = 1
        if isinstance(testCrossTalk,bool) and testCrossTalk:
            cTalkMat = np.random.rand(2 ** numQubits, 2 ** numQubits)
            for i in range(2 ** numQubits):
                cTalkMat[:, i] = cTalkMat[:, i] / (sum(cTalkMat[:, i] + 0 * np.random.random()))
            Tomo_Object.conf['Crosstalk'] = cTalkMat
        elif isinstance(testCrossTalk,np.ndarray):
            Tomo_Object.conf['Crosstalk'] = testCrossTalk
            testCrossTalk = True
        else:
            Tomo_Object.conf['Crosstalk'] = np.identity(2 ** numQubits)
        Tomo_Object.conf['UseDerivative'] = "TrUe"
        Tomo_Object.conf['Bellstate'] = testBell
        Tomo_Object.conf['DoErrorEstimation'] = errBounds
        if (testDrift):
            Tomo_Object.conf['DoDriftCorrection'] = "T"
        else:
            Tomo_Object.conf['DoDriftCorrection'] = "F"
        Tomo_Object.conf['Window'] = 1
        Tomo_Object.conf["Method"] = method
        if method.upper() == "HMLE":
            Tomo_Object.conf["Beta"] = 1/2
        if(TomographyNumber ==0):
            print("Running Test " + uniqueID(Tomo_Object),end=" .")

        tomo_input = Tomo_Object.getTomoInputTemplate()
        intensities = np.ones(tomo_input.shape[0])

        # set up measurements
        # measurements is an array of all the measurements for both classes mStates is only to help calculate Measurements
        if (Tomo_Object.getNumDetPerQubit() == 1):
            # input[:, np.arange(n_qubit+ 2, 3*n_qubit+ 2)]: measurements
            measurements_raw = tomo_input[:, np.arange(numQubits + 2, 3 * numQubits + 2)]
        else:
            # input[:, np.arange(2**n_qubit+ 2*n_qubit+ 1, 2**n_qubit+ 4*n_qubit+ 1)]: measurements
            measurements_raw = tomo_input[:, np.arange(2 ** numQubits + 2 * numQubits + 1,
                                              2 ** numQubits + 4 * numQubits + 1)]
        crosstalk = Tomo_Object.conf['Crosstalk']


        # NEW WAY
        measurements_densities = np.zeros((tomo_input.shape[0], Tomo_Object.getNumCoinc(), 2 ** numQubits, 2 ** numQubits), dtype=complex)
        for j in range(tomo_input.shape[0]):
            meas_basis_densities = np.zeros([2 ** numQubits, 2 ** numQubits, 2 ** numQubits]) + 0j
            meas_basis_pures = 1
            for k in range(numQubits):
                alpha = measurements_raw[j][2 * k]
                beta = measurements_raw[j][2 * k + 1]
                psi_transmit = np.array([alpha, beta])
                psi_reflect = np.array([np.conj(beta), np.conj(-alpha)])
                meas_pure = np.outer((np.array([1, 0])), psi_transmit) + np.outer((np.array([0, 1])), psi_reflect)
                # Changed from tensor_product to np.kron
                meas_basis_pures = np.kron(meas_basis_pures, meas_pure)
            for k in range(2 ** numQubits):
                meas_basis_densities[k, :, :] = np.outer(meas_basis_pures[:, k].conj().transpose(),
                                                         meas_basis_pures[:, k])
            for k in range(Tomo_Object.getNumCoinc()):
                for l in range(2 ** numQubits):
                    measurements_densities[j, k, :, :] = measurements_densities[j, k, :, :] + meas_basis_densities[l, :,
                                                                                              :] * crosstalk[k, l]
        # The old way of simulating the measurement
        # mStates = np.reshape(mStates, (mStates.shape[0], int(mStates.shape[1] / 2), 2))
        # if (test2Det):
        #     wavePlateArraysBasis = np.zeros((3, 2, 2), dtype=complex)
        #     wavePlateArraysBasis[0] = np.identity(2, dtype=complex)
        #     wavePlateArraysBasis[1] = np.array([[.7071, .7071], [.7071, -.7071]], dtype=complex)
        #     wavePlateArraysBasis[2] = np.array([[.7071, -.7071j], [.7071, .7071j]], dtype=complex)
        #     wavePlateArray = np.array([1], dtype=complex)
        # else:
        #     wavePlateArraysBasis = np.zeros((6, 2, 2), dtype=complex)
        #     wavePlateArraysBasis[0] = np.identity(2, dtype=complex)
        #     wavePlateArraysBasis[1] = np.array([[0, 1], [1, 0]], dtype=complex)
        #     wavePlateArraysBasis[2] = np.array([[.7071, .7071], [.7071, -.7071]], dtype=complex)
        #     wavePlateArraysBasis[3] = np.array([[.7071, -.7071], [.7071, .7071]], dtype=complex)
        #     wavePlateArraysBasis[4] = np.array([[.7071, -.7071j], [.7071, .7071j]], dtype=complex)
        #     wavePlateArraysBasis[5] = np.array([[.7071, .7071j], [.7071, -.7071j]], dtype=complex)
        #     wavePlateArray = np.array([1], dtype=complex)
        # for i in range(numQubits):
        #     wavePlateArray = np.kron(wavePlateArray, wavePlateArraysBasis)
        #
        # for i in range(len(mStates)):
        #     if (test2Det):
        #         for x in range(0, 2 ** (numQubits)):
        #             if (x == 1):
        #                 try:
        #                     mStates[i][1] = getOppositeState(mStates[i][1])
        #                 except:
        #                     mStates[i][0] = getOppositeState(mStates[i][0])
        #             elif (x == 2):
        #                 mStates[i][1] = getOppositeState(mStates[i][1])
        #                 mStates[i][0] = getOppositeState(mStates[i][0])
        #             elif (x == 3):
        #                 mStates[i][1] = getOppositeState(mStates[i][1])
        #             temp = mStates[i][0]
        #             for j in range(1, len(mStates[i])):
        #                 temp = np.kron(temp, mStates[i][j])
        #             measurements[i][x] = temp
        #         try:
        #             mStates[i][1] = getOppositeState(mStates[i][1])
        #         except:
        #             mStates[i][0] = getOppositeState(mStates[i][0])
        #         mStates[i][0] = getOppositeState(mStates[i][0])
        #     else:
        #         temp = mStates[i][0]
        #         for j in range(1, len(mStates[i])):
        #             temp = np.kron(temp, mStates[i][j])
        #         measurements[i] = temp
        if (TomographyNumber == 0):
            print(".", end="")

        ####################################
        # CREATE STATE AND SIM PROJECTIONS #
        ####################################

        # create random state
        if isinstance(randomStateDist,str):
            if randomStateDist == "density":
                startingRho = qLib.random_density_state(numQubits)
            elif randomStateDist == "pure":
                startingRho = qLib.random_pure_state(numQubits)
            elif randomStateDist == "bellstate":
                startingRho = qLib.random_bell_state(numQubits)
            else:
                raise ValueError('"'+ randomStateDist +'" is not a valid state distribution. '
                                'Possible values are ["density","pure","bellstate"]')
        elif isinstance(randomStateDist, np.ndarray):
                if abs(np.log2(randomStateDist.shape[0])-numQubits) < 10**-6:
                    if len(randomStateDist.shape) ==1:
                        startingRho = qLib.toDensity(randomStateDist)
                    else:
                        startingRho = randomStateDist
                else:
                    raise ValueError('"' + str(randomStateDist) + '" does not have the right dimensions for '+str(numQubits)+' Qubits.')
        else:
            raise ValueError('"' + str(randomStateDist) + '" is not a valid state distribution.')

        # counts
        numCounts = int(np.random.randint(np.ceil(5 * 2 ** numQubits), np.floor(50 * 2 ** numQubits)))

        # Testing setting
        for i in range(len(tomo_input)):
            if (test2Det):
                prob = np.zeros(2 ** numQubits, complex)
                for j in range(0, 2 * numQubits):
                    prob[j] = np.trace(np.matmul(measurements_densities[i,j], startingRho))
                    prob[j] = min(np.real(prob[j]), .99999999)
                prob = np.real(prob)
                tomo_input[i,2 * numQubits + 1: 2 ** numQubits + 2 * numQubits + 1] = np.random.multinomial(
                    numCounts, prob)
            else:
                prob = np.trace(np.matmul(measurements_densities[i,0], startingRho))
                prob = np.real(prob)
                tomo_input[i, numQubits + 1] = np.random.binomial(numCounts, min(prob, .99999999))

        if (testAccCorr and numQubits<3):
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
                acc[:, j] = np.prod(np.real(sings[:, tuple(index)]), axis=1) * (window[j] * 1e-9 / np.real(t)) ** (numQubits - 1)
            if (acc.shape != coinc.shape):
                acc = acc[:, 0]
            if (test2Det):
                tomo_input[:, np.arange(1, 2 * numQubits + 1)] = sings
                tomo_input[:,
                np.arange(2 * numQubits + 1, 2 ** numQubits + 2 * numQubits + 1)] = coinc + acc
            else:
                tomo_input[:, np.arange(1, numQubits + 1)] = sings
                tomo_input[:, numQubits + 1] = coinc + acc
            Tomo_Object.conf['Window'] = window
            tomo_input[:, 0] = t

        if (testDrift):
            intensities = np.random.random(intensities.shape)*2+.25
            if (test2Det):
                # tomo_input[:, np.arange(2*n_qubit+ 1, 2**n_qubit+ 2*n_qubit+ 1)]: coincidences
                coinc_range = np.arange(2 * numQubits + 1, 2 ** numQubits + 2 * numQubits + 1)
                for k in range(tomo_input.shape[0]):
                    tomo_input[k, coinc_range] = (intensities[k]* tomo_input[k, coinc_range]).real.astype(int)
            else:
                # tomo_input[:, n_qubit+ 1]: coincidences
                coinc_range = numQubits + 1
            for k in range(tomo_input.shape[0]):
                tomo_input[k, coinc_range] = (intensities[k] * tomo_input[k, coinc_range]).real.astype(int)

        if (TomographyNumber == 0):
            print(".")

        # Do tomography with settings
        try:
            start_time = time.time()
            myDensity, inten, myfVal = Tomo_Object.StateTomography_Matrix(tomo_input, intensities, method=method)
            end_time = time.time()

            Fidelity_with_Original = qLib.fidelity(startingRho, myDensity)
            Original_Purity = qLib.purity(startingRho)
            if (testBell):
                Tomo_Object.getBellSettings()
            if (Fidelity_with_Original < .8 and
                    not testCrossTalk and
                    np.average(Tomo_Object.getCoincidences()) > 10 and
                    method.upper() != "LINEAR"):
                print("-----------------------------")
                print("Low Fidelity of " + str(Fidelity_with_Original) + ". Avg counts per basis = " + str(
                    np.average(Tomo_Object.getCoincidences())))
        except:
            Original_Purity = -1
            Fidelity_with_Original = -1
            print("-----------------------------\n")
            print(traceback.format_exc() + "\n")
            numErrors += 1
            end_time = time.time()
            
        BigListOfTomographies.append([Tomo_Object, Fidelity_with_Original, Original_Purity, end_time - start_time])

    if (numErrors > 0):
        print("-----------------------------")
        print(str(numErrors)+' out of ' + str(nStates) + ' tomographys failed\n\n')
    else:
        print('No Issues!\n\n')
    return BigListOfTomographies


# Returns 1 if success
def saveRunsGeneral(numQubits, nStates,resultsFilePath="Results/results_GeneralData.csv",
                 randomStateDist="density", errBounds=0, testAccCorr=False, test2Det=False, testCrossTalk=False,
                 testBell=False, testDrift=False, saveData=True, method='MLE'):

    Tomographys = runTests(numQubits, nStates,randomStateDist=randomStateDist, errBounds=errBounds, testAccCorr=testAccCorr, 
                           test2Det=test2Det, testCrossTalk=testCrossTalk,testBell=testBell, testDrift=testDrift, method=method)

    numTotalErrors = 0

    if saveData:
        if os.path.isfile(resultsFilePath):
            printHeader=False
        else:
            printHeader=True

        # Open a file with access mode 'a'
        file_object = open(resultsFilePath, 'a')
    # Print details of each tomography to csv
    for t in Tomographys:
        [Tomo_Object, Fidelity_with_Original, Original_Purity, Total_Time] = t
        if Original_Purity <0:
            numTotalErrors +=1
        if saveData:
            # Set up dictionary
            dataRow = dict()
            dataRow['nqubits'] = Tomo_Object.conf["NQubits"]
            dataRow['method'] = Tomo_Object.conf["Method"]
            dataRow['fid_with_actual'] = Fidelity_with_Original
            dataRow['avg_coinc_per_meas'] = np.average(Tomo_Object.getCoincidences())
            dataRow['total_time'] = Total_Time
            dataRow['error_occurred'] = Original_Purity == -1
            dataRow['test_error'] = Tomo_Object.conf["DoErrorEstimation"]
            dataRow['test_acc'] = Tomo_Object.conf["DoAccidentalCorrection"]
            dataRow['test_det'] = Tomo_Object.conf["NDetectors"]
            dataRow['test_cross'] = np.any(Tomo_Object.conf["Crosstalk"] - np.eye(Tomo_Object.conf["Crosstalk"].shape[0]) > 1e-6)
            dataRow['test_bell'] = Tomo_Object.conf["Bellstate"]
            dataRow['test_drift'] = Tomo_Object.conf["DoDriftCorrection"]
            dataRow['purity'] = Original_Purity

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
    if saveData:
        # Close the file
        file_object.close()
    return numTotalErrors

def uniqueID(Tomo_Object):
    numQubits = Tomo_Object.conf["NQubits"]
    method = Tomo_Object.conf["Method"]
    errBounds = Tomo_Object.conf["DoErrorEstimation"]
    testAccCorr = Tomo_Object.conf["DoAccidentalCorrection"]
    numDet = Tomo_Object.conf["NDetectors"]

    if np.any(Tomo_Object.conf["Crosstalk"] - np.eye(Tomo_Object.conf["Crosstalk"].shape[0]) > 10 ** -6):
        testCrossTalk = True
    else:
        testCrossTalk = False

    testBell = Tomo_Object.conf["Bellstate"]
    testDrift = Tomo_Object.conf["DoDriftCorrection"]
    
    s = "N" + str(numQubits) + "-"
    s +="e" + str(errBounds) + "-"
    if (testAccCorr):
        s+="a1-"
    else:
        s += "a0-"
    if numDet == 2:
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
    s += "-" + method
    return s

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