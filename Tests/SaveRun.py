import numpy as np
import matplotlib.pyplot as plt
import QuantumTomography as qLib
import traceback
import time

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This class runs multiple tomographies with given settings. Create an instance of the SaveRun class 
and give it the args(or settings) to use. You can then use the run function and it will run  the tomographies"""
"""Data will not be deleted unless its by the user"""

"Attention! These tests run on the version that your environment uses. see readme for details"

class SaveRun():

    # todo:
    #   implement testAcc. Generate tomography with accidentals.

    counter = 0

    # private variables to store the data of each tomography
    allCounts = np.zeros(0,dtype=int)
    allFidels = np.zeros(0,dtype=float)
    allTimes = np.zeros(0,dtype=float)

    def __init__(self,args,Save_each_State=True,showErrors=True):

        [nBits, bounds, acc, det2, cross, bell, drift, nStates] = args

        """-------Settings---------"""
        self.numQubits = nBits # numeric >= 1
        self.nStates = nStates # numeric > 0
        self.test2Det = det2 # 0 or 1
        self.testCrossTalk = cross # 0 or 1
        self.errBounds = bounds # numeric >= 0
        self.testBell = bell # 0 or 1
        self.testDrift = drift # 0 or 1

        # not implemented
        self.testAccCorr = acc # 0 or 1

        # Turn this to false if you do not want an html output of each random test
        self.Save_each_State = Save_each_State

        # Turn this to false if you do not want show plot points where errors occured
        self.showErrors = showErrors

    # Handles cases where the settings are invalid. You can't do accidental correction with 1 qubit.
    def isValid(self):
        if(self.numQubits ==1 and (self.testAccCorr or self.testBell)):
            return False
        SaveRun.counter +=1
        return True

    # MAIN FUNCTION
    # ---------------
    # When you initialize a class you give it the settings. You then use this function to generate random
    # states, do tomography, and save the results
    def run(self):

        # csvSaver was a script to save the tomographies to a csv. It has been commented out
        # dataSaver = csvSaver(self.numQubits)

        # Declare Arrays to store the tomo data
        startingRhos = np.zeros((self.nStates, 2 ** self.numQubits, 2 ** self.numQubits), dtype=complex)
        myDensities = np.zeros((self.nStates, 2 ** self.numQubits, 2 ** self.numQubits), dtype=complex)
        myfVals = np.zeros((self.nStates), dtype=complex)
        myFidels = np.zeros((self.nStates))
        totalCounts = np.zeros((self.nStates))
        if (self.testCrossTalk):
            myCTalks = np.zeros((self.nStates, 2 ** self.numQubits, 2 ** self.numQubits), dtype=complex)

        # Check if the setting configuration is valid
        # Used to ignore settings that arent compatable
        success = False
        if not (self.isValid()):
            return 0

        # Decalare tomography class
        tomo = qLib.Tomography()

        # set up Conf Settings of the tomography class
        tomo.conf['NQubits'] = self.numQubits
        tomo.conf['Properties'] = ['concurrence', 'tangle', 'entanglement', 'entropy', 'linear_entropy', 'negativity']
        if(self.testAccCorr):
            tomo.conf['DoAccidentalCorrection'] = 1
        else:
            tomo.conf['DoAccidentalCorrection'] = 0
        if(self.test2Det):
            tomo.conf['NDetectors'] = 2
        else:
            tomo.conf['NDetectors'] = 1
        if(not self.testCrossTalk):
            tomo.conf['Crosstalk'] = np.identity(2 ** self.numQubits)
        tomo.conf['UseDerivative'] = 0
        tomo.conf['Bellstate'] = self.testBell
        tomo.conf['DoErrorEstimation'] = self.errBounds
        if(self.testDrift):
            tomo.conf['DoDriftCorrection'] = 1
        else:
            tomo.conf['DoDriftCorrection'] = 0
        tomo.conf['Window'] = 1
        tomo.conf['Efficiency'] = np.ones(2**self.numQubits)


        # Set up the tomo_input matrix
        tomo_input = tomo.getTomoInputTemplate()

        # Fill the tomo_input measurement columns
        # The measurments are based on using wave plates.
        #   WavePlate basis is a list of wavePlate matricies.
        #   WavePlateArray is a list of the overall gate that is applied to the state.
        #   The model is the random state goes through the wavePlateArray[i] gate, and then it is
        #   projected onto 1. This is the same as the state being projected onto measurments[i]
        tomo_input = tomo.getTomoInputTemplate()
        if(self.test2Det):
            measurements = np.zeros((len(tomo_input),2**self.numQubits,2**(self.numQubits)),dtype=complex)
        else:
            measurements = np.zeros((len(tomo_input),2**self.numQubits),dtype=complex)
        if (tomo.getNumDetPerQubit() == 1):
            # input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
            mStates = tomo_input[:,np.arange(self.numQubits + 2, 3 * self.numQubits + 2)]
        else:
            # input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements
            mStates = tomo_input[:,np.arange(2 ** self.numQubits + 2 * self.numQubits + 1, 2 ** self.numQubits + 4 * self.numQubits + 1)]
        mStates = np.reshape(mStates, (mStates.shape[0],int(mStates.shape[1]/2), 2))
        if(self.test2Det):
            wavePlateArraysBasis = np.zeros((3, 2, 2), dtype=complex)
            wavePlateArraysBasis[0] = np.identity(2, dtype=complex)
            wavePlateArraysBasis[1] = np.array([[.7071, .7071], [.7071, -.7071]], dtype=complex)
            wavePlateArraysBasis[2] = np.array([[.7071, -.7071j], [.7071, .7071j]], dtype=complex)
            wavePlateArray = np.array([1], dtype=complex)
        else:
            wavePlateArraysBasis = np.zeros((6,2,2),dtype=complex)
            wavePlateArraysBasis[0] = np.identity(2,dtype=complex)
            wavePlateArraysBasis[1] = np.array([[0,1],[1,0]],dtype=complex)
            wavePlateArraysBasis[2] = np.array([[.7071,.7071],[.7071,-.7071]],dtype=complex)
            wavePlateArraysBasis[3] = np.array([[.7071,-.7071],[.7071,.7071]],dtype=complex)
            wavePlateArraysBasis[4] = np.array([[.7071,-.7071j],[.7071,.7071j]],dtype=complex)
            wavePlateArraysBasis[5] = np.array([[.7071,.7071j],[.7071,-.7071j]],dtype=complex)
            wavePlateArray = np.array([1],dtype=complex)
        for i in range(self.numQubits):
            wavePlateArray = np.kron(wavePlateArray,wavePlateArraysBasis)


        # set up mState array. Contains all the measurements
        #   mState[i] is a list of measurements on each qubit.
        #   mState[i,0] is measurment on the first qubit
        #   mState[i,1] is measurment on the second qubit
        # also set up measurments array. List of overall measurent
        #   measuremets[i] is the total quantum state the state is measured on. It
        #   is the kronecker product of mState[i,0] x ... x mState[i,j] x ... x mState[i,n]
        for i in range(len(mStates)):
            if(self.test2Det):
                for x in range(0,2**(self.numQubits)):
                    if(x == 1):
                        try:
                            mStates[i][1] = getOppositeState(mStates[i][1])
                        except:
                            mStates[i][0] = getOppositeState(mStates[i][0])
                    elif(x == 2):
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
                for j in range(1,len(mStates[i])):
                    temp = np.kron(temp, mStates[i][j])
                measurements[i] = temp

                # measurements is an array of all the measurements for both classes mStates is only to help calculate Measurements


        '''###########'''
        ''' Main Loop '''
        '''###########'''
        # This loop creates random states and runs tomography. It then stores the results in
        # the matrices defined above.

        # keep a running track of number of errors
        numErrors = 0
        #create states and do tomo
        for x in range(self.nStates):
            # if crosstalk is enabled then create a random crosstalk matrix
            cTalkMat = 0
            if(self.testCrossTalk):
                cTalkMat = np.random.rand(2**self.numQubits,2**self.numQubits)
                for i in range(2**self.numQubits):
                    cTalkMat[:,i] = cTalkMat[:,i]/(sum(cTalkMat[:,i]+0*np.random.random()))
                tomo.conf['Crosstalk'] = cTalkMat
                myCTalks[x] = cTalkMat

            # counts
            numCounts = np.zeros(tomo_input.shape[0],dtype=int)
            exponentForCounts = np.random.randint(1, 6)
            for c in range(tomo_input.shape[0]):
                numCounts[c] = (1+np.random.random())*10**exponentForCounts

            # create random state
            state = qLib.random_pure_state(self.numQubits)
            startingRho = qLib.toDensity(state)

            # Run the state through the waveplates and then project onto h basis which is the |0^n> state.
            #   This loop genrates the counts probabalistically
            for i in range(len(tomo_input)):
                # state goes through wave plates
                newState = np.matmul(wavePlateArray[i], state)
                # state goes through beam splitter and we measure the H counts
                if (self.testCrossTalk):
                    newState = np.matmul(cTalkMat, newState)
                hBasis = np.zeros(2 ** self.numQubits, dtype=complex)
                hBasis[0] = 1
                if (self.test2Det):
                    prob = np.zeros(2 ** self.numQubits, complex)
                    for j in range(0, 2 * self.numQubits):
                        hBasis = np.zeros(2**self.numQubits, dtype=complex)
                        hBasis[j] = 1
                        prob[j] = np.dot(hBasis, newState)
                        prob[j] = min(prob[j] * prob[j].conj(), .99999999)
                    prob = np.array(prob, dtype=float)
                    tomo_input[i, 2 * self.numQubits + 1: 2 ** self.numQubits + 2 * self.numQubits + 1] = np.random.multinomial(numCounts[i], prob)
                else:
                    prob = np.dot(hBasis, newState)
                    prob = prob * prob.conj()
                    tomo_input[i, self.numQubits + 1] = np.random.binomial(numCounts[i],min(prob,.99999999))
                # input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences

            # Generate Data for accidental corrections.
            # TODO: double check this is correct
            if (self.testAccCorr):
                # acc[:, j] = np.prod(np.real(sings[:, idx]), axis=1) * (window[j] * 1e-9 / np.real(t)) ** (nbits - 1)
                if(self.test2Det):
                    sings = tomo_input[:, np.arange(1, 2*self.numQubits + 1)]
                    coinc = tomo_input[:,np.arange(2 * self.numQubits + 1, 2 ** self.numQubits + 2 * self.numQubits + 1)]
                    window = np.ones(4)*100000
                    n_coinc = 2**self.numQubits
                else:
                    sings = tomo_input[:, np.arange(1, self.numQubits + 1)]
                    coinc = tomo_input[:, self.numQubits+1]
                    window = [100000]
                    n_coinc = 1
                sings = 2000*np.random.random(sings.shape)
                t = np.random.random(tomo_input.shape[0])+1
                acc = np.zeros_like(coinc)
                if (len(acc.shape) == 1):
                    acc = acc[:, np.newaxis]

                scalerIndex = np.concatenate((np.ones(self.numQubits - 2), [2, 2]))
                additiveIndex = np.array([0, 1])
                for j in range(2, self.numQubits):
                    additiveIndex = np.concatenate(([2 * j], additiveIndex))
                for j in range(n_coinc):
                    index = bin(j).split("b")[1]
                    index = "0" * (self.numQubits - len(index)) + index
                    index = [int(char) for char in index]
                    index = index * scalerIndex + additiveIndex
                    index = np.array(index, dtype=int)
                    acc[:, j] = np.prod(np.real(sings[:, tuple(index)]), axis=1) * (window[j] * 1e-9 / np.real(t)) ** (self.numQubits - 1)
                if (acc.shape != coinc.shape):
                    acc = acc[:, 0]
                if (self.test2Det):
                    tomo_input[:, np.arange(1, 2 * self.numQubits + 1)] = sings
                    tomo_input[:,np.arange(2 * self.numQubits + 1, 2 ** self.numQubits + 2 * self.numQubits + 1)] = coinc + acc
                else:
                    tomo_input[:, np.arange(1, self.numQubits + 1)] = sings
                    tomo_input[:, self.numQubits + 1] = coinc + acc
                tomo.conf['Window'] = window
                tomo_input[:, 0] = t

            # Generate Data for drift corrections.
            if (self.testDrift):
                intensity = numCounts
                if (self.test2Det):
                    # tomo_input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
                    for k in range(tomo_input.shape[0]):
                        tomo_input[k, np.arange(2 * self.numQubits + 1, 2 ** self.numQubits + 2 * self.numQubits + 1)] = int((intensity[k])*tomo_input[k, np.arange(2*self.numQubits+1, 2**self.numQubits+2*self.numQubits+1)])
                else:
                    # tomo_input[:, n_qubit+1]: coincidences
                    for k in range(tomo_input.shape[0]):
                        tomo_input[k, self.numQubits+1] = int((intensity[k])*tomo_input[k, self.numQubits+1])


            # State_Tomography CALL
            # This is calls the function state_Tomography and times how long it took.
            # It can also handle cases where the tomography failed and encountered an error
            errorMessage = ""

            try:
                start_time = time.time()
                myDensitie, inten, myfVal = tomo.state_tomography(tomo_input, numCounts)
                tomographyTime = (time.time() - start_time)
                myFidel = qLib.fidelity(startingRho,myDensitie)
                if(self.testBell):
                    tomo.getBellSettings(myDensitie)
                    tomo.getProperties(myDensitie)
            except:
                myDensitie = [[0]]
                inten = 0
                myfVal = -1
                myFidel = -1
                errorMessage = "<pre>\n"+traceback.format_exc()+"\n</pre>"



            # Save data from the last tomography and put it into our data matrices
            # dataSaver.addData(sum(numCounts),myFidel,tomographyTime)
            startingRhos[x] = startingRho
            myDensities[x] = myDensitie
            myfVals[x] = myfVal
            myFidels[x] = myFidel
            totalCounts[x] = sum(numCounts)

            # We will say we encountered an error if our fidelity is less than .8
            if(myFidel < .8):
                numErrors += 1

        '''##############'''
        ''' Save Results '''
        '''##############'''
        # This code generates an html and pdf.
        # dataSaver.saveData() was previously used but has been commented out

        # Create Graph of Fidelities
        fig = plt.figure()
        fig.clf()

        if not self.showErrors:
            # Change zeros to nan so plt.plot won't plot them
            numNANS = 0
            for e in range(len(myFidels)):
                if myFidels[e] == -1:
                    myFidels[e] = float('nan')
                    numNANS +=1
            plt.text(-.04, -.085, "Showing " + str(len(myFidels) - numNANS) + " out of " + str(len(myFidels)) + " trials. Num Errors: " + str(numNANS))
            plt.subplots_adjust(bottom=0.25)

        plt.plot(np.log10(totalCounts), myFidels, '.b')
        plt.title('Fidelities')

        plt.xlabel("Log(Counts) base 10")
        plt.ylabel("Fidelity")
        plt.savefig("Results/Fidel_Graph_"+self.uniqueID()+".png")


        ''' CREATE HTML PAGE '''
        # FOROREPLACE is a big string with most of the html code. Some stuff like the header is
        # declared in the variable fff. The following code goes through the results and adds
        # to this string. At the end the string is written to html file in the results folder


        # This chunk shows what settings were used in the top right.
        FORREPLACE = '<table style="margin: auto;width: 70%;">'
        FORREPLACE += '<tr><th>Settings used</th><th></th></tr>'
        FORREPLACE += '<tr><td><ul>'
        # Test ID
        FORREPLACE += '<li>'
        FORREPLACE += '<b>Test ID : </b>'
        FORREPLACE += self.uniqueID()
        FORREPLACE += '</li>'
        # num qubits
        FORREPLACE += '<li>'
        FORREPLACE += '<b>Number of Qubits : </b>'
        FORREPLACE += str(self.nStates)
        FORREPLACE += '</li>'
        # num States
        FORREPLACE += '<li>'
        FORREPLACE += '<b>Number of States : </b>'
        FORREPLACE += str(self.numQubits)
        FORREPLACE += '</li>'
        #err correction
        FORREPLACE += '<li>'
        FORREPLACE += '<b>Number of errCorr : </b>'
        FORREPLACE += str(self.errBounds)
        FORREPLACE += '</li>'
        if (self.test2Det):
            FORREPLACE += '<li>'
            FORREPLACE += '<b>test2Det : </b> True'
            FORREPLACE += '</li>'
        if (self.testCrossTalk):
            FORREPLACE += '<li>'
            FORREPLACE += '<b>testCrossTalk : </b> True'
            FORREPLACE += '</li>'
        if (self.testBell):
            FORREPLACE += '<li>'
            FORREPLACE += '<b>testBell : </b> True'
            FORREPLACE += '</li>'
        if (self.testAccCorr):
            FORREPLACE += '<li>'
            FORREPLACE += '<b>testAccCorr : </b> True'
            FORREPLACE += '</li>'
        if (self.testDrift):
            FORREPLACE += '<li>'
            FORREPLACE += '<b>testDrift : </b> True'
            FORREPLACE += '</li>'

        # adds the fidelity graph
        FORREPLACE += '</ul></td><td><img src="Fidel_Graph_'+self.uniqueID()+'.png" style="width:80%;height:auto;"></td></tr></table>'

        if(self.Save_each_State):
            for j in range(self.nStates):
                '''Real State'''
                FORREPLACE += '<table class="data">'
                FORREPLACE += '<tr><th colspan="2">Actual Densities</th><th>Details</th>'

                FORREPLACE += '<th colspan="2">Caclulated Densities</th>'
                if (self.testCrossTalk):
                    FORREPLACE += '<th colspan="2">CrossTalk</th>'
                FORREPLACE += '</tr>'
                FORREPLACE += '<tr>'
                FORREPLACE += '<td colspan="2">'
                FORREPLACE += qLib.matrixToHTML(startingRhos[j])
                FORREPLACE += '</td>'
                
                '''My tomography'''
                FORREPLACE += '<td>'
                FORREPLACE += printfValsFidelity(myfVals[j], myFidels[j],totalCounts[j])
                FORREPLACE += '</td>'
                FORREPLACE += '<td colspan="2">'
                FORREPLACE += qLib.matrixToHTML(myDensities[j])
                FORREPLACE += '</td>'

                # Print crosstalk matrix
                if (self.testCrossTalk):
                    FORREPLACE += '<td colspan="2">'
                    FORREPLACE += qLib.matrixToHTML(myCTalks[j])
                    FORREPLACE += '</td>'

                FORREPLACE += '</tr>'
            FORREPLACE += '</table>'

        # Edit and Save HTML file
        fff = '<html><head><title>Data</title><link href="styleSheet.css" rel="stylesheet" type="text/css"></head><body>TOREPLACE</body></html>'
        fff = fff.replace('TOREPLACE', str(FORREPLACE))
        with open('Results/RandomDataOutPut_'+self.uniqueID()+'.html', 'w') as ff:
            ff.write(fff)
            ff.close()


        # Print out details of what was tested in the console
        print("-----------------------")
        print("Test  "+str(SaveRun.counter)+": "+self.uniqueID())
        print("Number of Qubits Bits: "+ str(self.numQubits))
        print("Number of errCorr: " + str(self.errBounds))
        if (self.test2Det):
            print("test2Det:True")
        else:
            print("test2Det:False")
        if (self.testCrossTalk):
            print("testCrossTalk:True")
        else:
            print("testCrossTalk:False")
        if (self.testBell):
            print("testBell:True")
        else:
            print("testBell:False")
        if (self.testAccCorr):
            print("testAccCorr:True")
        else:
            print("testAccCorr:False")
        if (self.testDrift):
            print("testDrift:True")
        else:
            print("testDrift:False")

        print("Result : COMPLETED " + "[" + str(numErrors) + "/" + str(self.nStates) + "] Errors")
        print("-----------------------\n")

        # THIS IS THE END OF THE RUN FUNCTION
        return 1


    # The following are helper functions used during run. Some of these may not be used.


    # Returns a string based on the current args
    def uniqueID(self):
        s = "N"+str(self.numQubits) +"-"
        s +="e"+str(self.errBounds) +"-"
        if (self.testAccCorr):
            s+="a1-"
        else:
            s += "a0-"
        if (self.test2Det):
            s+="d1-"
        else:
            s += "d0-"
        if (self.testCrossTalk):
            s+="c1-"
        else:
            s += "c0-"
        if (self.testBell):
            s+="b1-"
        else:
            s += "b0-"
        if (self.testDrift):
            s+="dr1"
        else:
            s += "dr0"
        return s
#
# def selectionSort(results):
#     # Traverse through all array elements
#     sortedResults = results.copy()
#     for i in range(0, len(results)):
#         # Find the minimum element in remaining unsorted array
#         min_idx = i
#         for j in range(i + 1, len(sortedResults)):
#             if sortedResults[min_idx] > sortedResults[j]:
#                 min_idx = j
#         # Swap the found minimum element with the first element
#         if (i != min_idx):
#             sortedResults[min_idx], sortedResults[i] = swap(sortedResults[min_idx], sortedResults[i])
#             for x in sortedResults:
#                 x[min_idx], x[i] = swap(x[min_idx], x[i])
#         return sortedResults

def printfValsFidelity(fVal,fid,Counts):
    lowerBoundOfFidelity = .8
    colorR = mapVal(fid,1,lowerBoundOfFidelity,0,255) #red
    res = '<ul class="Details" style="color:rgb(' + str(colorR) + ',0,0);">'
    if (fVal == -1):
        res += "<li><b>ERROR OCCURED<b></li>"
    else:
        # Counts
        res += '<li>'
        res += '<b>Counts : </b>' + str(int(Counts))
        res += '</li>'

        # fVal
        res += '<li>'
        res += '<b>fVal : </b>' + str(round(np.real(fVal), 6))
        res += '</li>'

        # Fidelity
        res += '<li>'
        res += '<b>Fidelity : </b>' + str(round(fid, 6))
        res += '</li>'

    res += str('</ul>')
    return res

#
# def swap(ele1, ele2):
#     temp = ele1.copy()
#     ele1 = ele2.copy()
#     ele2 = temp.copy()
#     return ele1, ele2

# hard coded mapping of states. Get the other state in the specified pair
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

# class csvSaver():
#     def __init__(self,numQubits):
#         self.n = numQubits
#         self.previousData = self.loadData()
#         self.count = 0
#
#         # determines how many states it should hold before saving it to the csv file.
#         saveNumber = ["Error",500,100,5]
#         self.nStates = 1
#         try:
#             self.nStates = saveNumber[self.n]
#         except :
#             self.nStates = 1
#         self.newData = np.zeros((3,self.nStates))
#
#     def loadData(self):
#         try:
#             return np.loadtxt('Results/myLibrary'+str(self.n)+'.csv', delimiter=',')
#         except:
#             return np.zeros((3,0))
#     def saveData(self):
#         finalData = np.concatenate((self.previousData,self.newData[:, ~np.all(self.newData == 0, axis=0)]), axis=1)
#         np.savetxt('Results/myLibrary'+str(self.n)+'.csv',finalData, delimiter=',')
#         self.previousData = finalData
#         self.newData = np.zeros((3, self.nStates))
#     def addData(self,count,fidel,time):
#         if(fidel!=-1):
#             self.newData[:,self.count] = [count,fidel,time]
#             self.count +=1
#             if(self.count>= self.nStates):
#                 self.saveData()
#                 self.count = 0
def mapVal(value,x1,x2,y1,y2,f = 1):
    newVal = value-x1

    newVal = (y2-y1)/(x2-x1) * newVal
    newVal += y1
    return newVal