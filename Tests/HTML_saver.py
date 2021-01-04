import numpy as np
import scipy.io
import sys
import matplotlib.pyplot as plt
sys.path.append("../Quantum-Tomography_Git")
import QuantumTomography as qLib

numQubits = 1
nStates = 30
testWebsite = False
"""----------------"""
testAccCorr = False
test2Det = False
testCrossTalk = True
testDrift = False

# Does not perform any tests.
testPerfectTomo = True
# Performes tests but does not use them in the calculations
testNoCorr = True
testBell = False
errBounds = 0


settingsArray = [testAccCorr,test2Det,testCrossTalk,testDrift,testWebsite]
settingsArrayName = ['testAccCorr','test2Det','testCrossTalk','testDrift','testWebsite']

tomo = qLib.Tomography()

'''pReal    totCounts       fVals       pCalculated     OptimizationNumber  meCalculated'''
'''eigStates                fS          eigStates       Fidelity            eigStates '''
''''''

def findM(meas):
    for x in range(0,len(data['projOps'])):
        if((np.linalg.norm(data['projOps'][x] - meas)<.0001)):
            return x;
    return -1
def doTomography(counts):
    if( int(np.log2(data['pReal'].shape[3]) ==2)):
        for x in range(0,36):
            measurement = np.kron(tomo_input[x,4:6],tomo_input[x,6:8])
            measurement = measurement[np.newaxis].conj().T @ measurement[np.newaxis]
            index = findM(measurement)
            tomo_input[x,3] = counts[index]
    else:
        tomo_input[0, 2] = counts[0]
        tomo_input[1, 2] = counts[1]
        tomo_input[2, 2] = counts[2]
        tomo_input[3, 2] = counts[3]
        tomo_input[4, 2] = counts[5]
        tomo_input[5, 2] = counts[4]

    if(testWebsite):
        [rho, intens, fval] = main_web_tester.calcTomography(tomo_input, intensity)
    else:
        [rho, intens, fval] = tomo.state_tomography(tomo_input, intensity)
    # rho = rho.conjugate()
    return [rho, intens, fval]

#Maps a range of numbers to another range of numbers
def map(value,x1,x2,y1,y2,f = 1):
    newVal = value-x1

    newVal = (y2-y1)/(x2-x1) * newVal
    newVal += y1
    return newVal

def printMatrix(M):
    s = np.shape(M)
    res = str('<table class="densityMat">')
    for i in range(s[0]):
        res += ' <tr>'
        for j in range(s[1]):
            res += '<td>' + str(round(M[i,j].real, 6)) + '<div class="purp">+</div><BR>' + str(round(M[i,j].imag, 6))
            res += '<div class="purp">j</div></td>'
        res += '</tr>'
    res += '</table>'
    d, v = np.linalg.eig(M)
    colorR = 0
    if(any(d < 0) or any(d > 1)):
        colorR = 255
    eigenVals = '<h5 class="eigVals" style="color:rgb(' + str(colorR) + ',0,0);">Eigen Values: </h5><p class="eigVals">'
    for x in range(0,len(d)):
        eigenVals += str(round(d[x].real, 5))
        if(abs(d[x].imag)>.00001):
            eigenVals += '<div class="purp">+</div>'
            eigenVals += str(round(d[x].imag, 5))
            eigenVals += '<div class="purp">j</div>'
        eigenVals +=" , "
    eigenVals = str(eigenVals)[0:len(str(eigenVals))-2]
    eigenVals += "<p>"
    res += eigenVals
    return res
def printfValsFidelity(fVal,fid,fs = -1):
    colorR = map(fid,1,lowerBoundOfFidelity,0,255) #red
    res = '<ul class="fVals">'
    res += '<li>'
    if(fVal == -1):
        res+= "ERROR OCCURED"
    else:
        res += '<div class="title">fVal: </div>'
        res += str(round(fVal, 6))
    res += '</li>'
    res += '<li style="color:rgb(' + str(colorR) + ',0,0);">'
    res += '<div class="title">Fidelity: </div><div class="inline">'
    res += str(round(fid, 6))
    res += '</div></li>'
    if(fs !=-1):
        colorR = map(fs,1,lowerBoundOfFidelity,0,255) #red
        res += '<li style="color:rgb(' + str(colorR) + ',0,0);">'
        res += '<div class="title">Chris\'s Fid Calc: </div><div class="inline">'
        res += str(round(fs, 6))
        res += '</div></li>'
    res += str('</ul>')
    return res

def swap(ele1,ele2):
    temp = ele1.copy()
    ele1 = ele2.copy()
    ele2 = temp.copy()
    return ele1,ele2

def indexOf(array, item):
    for x in range(0,array.shape[0]):
        if all((np.round(array[x],7) == np.round(item,7)).flatten()):
            return x
    raise Exception('could not find element')

def selectionSort(fids,array):
    # Traverse through all array elements
    for i in range(0,len(fids)):
        # Find the minimum element in remaining
        # unsorted array
        min_idx = i
        for j in range(i + 1, len(fids)):
            if fids[min_idx] > fids[j]:
                min_idx = j
        # Swap the found minimum element with the first element
        if(i!= min_idx):
            fids[min_idx], fids[i] = swap(fids[min_idx], fids[i])
            for x in array:
                x[min_idx], x[i] = swap(x[min_idx], x[i])
# oldData = [myFidels.copy(), startingRhos.copy(), totalCounts.copy(), myDensities.copy(), myfVals.copy()]


def checkSort(newData,oldData):
    for x in range(0,newData[1].shape[0]):
        index = indexOf(oldData[1],newData[1][x])
        # if( all((np.round(oldData[0][index],7) != np.round(newData[0][x],7)).flatten()) or \
        #     all((np.round(oldData[2][index],7) != np.round(newData[2][x],7)).flatten()) or \
        #     all((np.round(oldData[3][index],7) != np.round(newData[3][x],7)).flatten()) or \
        #     all((np.round(oldData[4][index],7) != np.round(newData[4][x],7)).flatten())):
        if (all((oldData[0][index] != newData[0][x]).flatten()) or \
                all((oldData[2][index] != newData[2][x]).flatten()) or \
                all((oldData[3][index] != newData[3][x]).flatten()) or \
                all((oldData[4][index] != newData[4][x]).flatten())):
            raise Exception('Sorting did not Maintain order')

def getOppositeState(psi):
    #Horizontal
    if(all(psi == np.array([1,0],dtype=complex))):
        return np.array([0,1],dtype=complex)
    if (all(psi == np.array([0, 1],dtype=complex))):
        return np.array([1, 0],dtype=complex)
    #Diagional
    if (all(psi == np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2))],dtype=complex))):
        return np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2))],dtype=complex)
    if (all(psi == np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2))],dtype=complex))):
        return np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2))],dtype=complex)
    #Circle
    if (all(psi == np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j],dtype=complex))):
        return np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j],dtype=complex)
    if (all(psi == np.array([(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j],dtype=complex))):
        return np.array([(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j],dtype=complex)
    else:
        raise Exception('State Not Found getOppositeState')
lowerBoundOfFidelity = .9

# what arrays we need
# startingRhos,myfVals,myFidels,myDensities,totalCounts
startingRhos = np.zeros((nStates, 2 ** numQubits, 2 ** numQubits), dtype=complex)
myDensities = np.zeros((nStates, 2 ** numQubits, 2 ** numQubits), dtype=complex)
myfVals = np.zeros((nStates), dtype=complex)
myFidels = np.zeros((nStates))

totalCounts = np.zeros((nStates))
AtotalCounts = np.zeros((nStates))

if(testCrossTalk):
    myCTalks = np.zeros((nStates, 2 ** numQubits, 2 ** numQubits), dtype=complex)

#set up settings for no Correction data
if(testPerfectTomo):
    perfectDensities = np.zeros((nStates, 2 ** numQubits, 2 ** numQubits), dtype=complex)
    perfectFVals = np.zeros((nStates), dtype=complex)
    perfectFidels = np.zeros((nStates))

    # set up settings for perfect tomo class
    tomoPerfect = qLib.Tomography()
    intensityPerfect = np.ones(6 ** numQubits)
    tomoPerfect.conf['NQubits'] = numQubits
    tomoPerfect.conf['Properties'] = ['concurrence', 'tangle', 'entanglement', 'entropy', 'linear_entropy', 'negativity']
    tomoPerfect.conf['DoAccidentalCorrection'] = 0
    if (test2Det):
        tomoPerfect.conf['NDetectors'] = 2
    else:
        tomoPerfect.conf['NDetectors'] = 1
    tomoPerfect.conf['Crosstalk'] = np.identity(2 ** numQubits)
    tomoPerfect.conf['UseDerivative'] = 0
    tomoPerfect.conf['Bellstate'] = 0
    tomoPerfect.conf['DoErrorEstimation'] = 0
    tomoPerfect.conf['DoDriftCorrection'] = 0
    tomoPerfect.conf['Window'] = 0
    tomoPerfect.conf['Efficiency'] = np.ones(2**numQubits)
    tomo_inputPerfect = tomoPerfect.getTomoInputTemplate()

#set up settings for default tomo class
intensity = np.ones(6 ** numQubits)
tomo.conf['NQubits'] = numQubits
tomo.conf['Properties'] = ['concurrence', 'tangle', 'entanglement', 'entropy', 'linear_entropy', 'negativity']
if(testAccCorr):
    tomo.conf['DoAccidentalCorrection'] = 1
else:
    tomo.conf['DoAccidentalCorrection'] = 0
if(test2Det):
    tomo.conf['NDetectors'] = 2
else:
    tomo.conf['NDetectors'] = 1
if(not testCrossTalk):
    tomo.conf['Crosstalk'] = np.identity(2 ** numQubits)
tomo.conf['UseDerivative'] = 0
tomo.conf['Bellstate'] = testBell
tomo.conf['DoErrorEstimation'] = errBounds
if(testDrift):
    tomo.conf['DoDriftCorrection'] = 1
else:
    tomo.conf['DoDriftCorrection'] = 0
tomo.conf['Window'] = 1
tomo.conf['Efficiency'] = np.ones(2**numQubits)

tomo_input = tomo.getTomoInputTemplate()

#set up measurements
if(test2Det):
    measurements = np.zeros((len(tomo_input),2**numQubits,2**(numQubits)),dtype=complex)
else:
    measurements = np.zeros((len(tomo_input),2**numQubits),dtype=complex)
if (tomo.getNumDetPerQubit() == 1):
    # input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
    mStates = tomo_input[:,np.arange(numQubits + 2, 3 * numQubits + 2)]
else:
    # input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements
    mStates = tomo_input[:,np.arange(2 ** numQubits + 2 * numQubits + 1, 2 ** numQubits + 4 * numQubits + 1)]
mStates = np.reshape(mStates, (mStates.shape[0],int(mStates.shape[1]/2), 2))
if(test2Det):
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
for i in range(numQubits):
    wavePlateArray = np.kron(wavePlateArray,wavePlateArraysBasis)




# I should put some of this code into the to Density function so that to puts multiple
for i in range(len(mStates)):
    if(test2Det):
        for x in range(0,2**(numQubits)):
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
''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
#create states and do tomo
for x in range(nStates):

    #to implement
    # testDrift
    # testAcc

    #set up crossyalk
    if(testCrossTalk):
        cTalkMat = np.random.rand(2**numQubits,2**numQubits)
        for i in range(2**numQubits):
            cTalkMat[:,i] = cTalkMat[:,i]/(sum(cTalkMat[:,i]+0*np.random.random()))

        tomo.conf['Crosstalk'] = cTalkMat
        myCTalks[x] = cTalkMat
    numCounts = int(np.random.random()*10**np.random.randint(4,5))

    # create random state
    state = np.random.beta(.5,.5,2**numQubits)+1j*np.random.beta(.5,.5,2**numQubits)
    state = state / np.sqrt(np.dot(state,state.conj()))
    startingRhos[x] = qLib.toDensity(state)

    # testing perfect
    if(testPerfectTomo):
        for i in range(len(tomo_input)):
            # state goes through wave plates
            newState = wavePlateArray[i] @ state
            # state goes throug beam splitter and we measure the H counts
            hBasis = np.zeros(2**numQubits, dtype=complex)
            if (test2Det):
                prob = np.zeros(2**numQubits,complex)
                for j in range(0,2*numQubits):
                    hBasis = np.zeros(2**numQubits, dtype=complex)
                    hBasis[j] = 1
                    prob[j] = np.dot(hBasis, newState)
                    prob[j] = min(prob[j] * prob[j].conj(),.99999999)
                prob = np.array(prob,dtype=float)
                tomo_inputPerfect[i, 2*numQubits+1: 2**numQubits+2*numQubits+1] = np.random.multinomial(numCounts, prob)
            else:
                hBasis[0] = 1
                prob = np.dot(hBasis, newState)
                prob = prob * prob.conj()
                tomo_inputPerfect[i, numQubits + 1] = np.random.binomial(numCounts, min(prob, .99999999))
        # Run perfect tomo
        try:
            [perfectDensities[x], totalCounts[x], perfectFVals[x]] = tomoPerfect.state_tomography(tomo_inputPerfect, intensity)
            perfectFidels[x] = fidelity(startingRhos[x],perfectDensities[x])
        except:
            perfectFVals[x] = -1
            perfectFidels[x] = 0

    #Testing settinga
    for i in range(len(tomo_input)):
        # state goes through wave plates
        newState = wavePlateArray[i] @ state
        # state goes through beam splitter and we measure the H counts
        if (testCrossTalk):
            newState = cTalkMat @ newState
        hBasis = np.zeros(2**numQubits, dtype=complex)
        if(test2Det):
            prob = np.zeros(2 ** numQubits, complex)
            for j in range(0, 2 * numQubits):
                hBasis = np.zeros(numQubits**2, dtype=complex)
                hBasis[j] = 1
                prob[j] = np.dot(hBasis, newState)
                prob[j] = min(prob[j] * prob[j].conj(), .99999999)
            prob = np.array(prob, dtype=float)
            tomo_input[i, 2 * numQubits + 1: 2 ** numQubits + 2 * numQubits + 1] = np.random.multinomial(
                numCounts, prob)
        else:
            hBasis[0] = 1
            prob = np.dot(hBasis, newState)
            prob = prob * prob.conj()
            tomo_input[i, numQubits + 1] = np.random.binomial(numCounts,min(prob,.99999999))
        # input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
    if (testAccCorr):
        if(test2Det):
            # tomo_input[:, np.arange(1, 2 * n_qubit + 1)]: singles
            pass
        else:
            # tomo_input[:, np.arange(1, n_qubit + 1)]: singles
            pass
    if (testDrift):
        if (test2Det):
            # tomo_input[:, 0]: times
            pass
        else:
            # tomo_input[:, 0]: times
            pass

    # Do tomography with settings
    try:
        [myDensities[x], totalCounts[x], myfVals[x]] = tomo.state_tomography(tomo_input, intensity)

        # what arrays we need
        # startingRhos,myfVals,myFidels,myDensities,totalCounts

        AtotalCounts[x] = numCounts
        myFidels[x] = fidelity(startingRhos[x],myDensities[x])
        tomo.getProperties(myDensities[x])
        if(testBell):
            tomo.getBellSettings(myDensities[x])
    except:
        AtotalCounts[x] = numCounts
        myfVals[x] = -1
        myFidels[x] = 0

# Sort data by fidelity
if(testCrossTalk):
    if(testPerfectTomo):
            selectionSort(myFidels, [startingRhos, totalCounts, myDensities, myfVals,perfectFidels,perfectDensities,perfectFVals,myCTalks])
    else:
            selectionSort(myFidels, [startingRhos, totalCounts, myDensities, myfVals,myCTalks])
else:
    if (testPerfectTomo):
            selectionSort(myFidels, [startingRhos, totalCounts, myDensities, myfVals, perfectFidels, perfectDensities,
                                     perfectFVals])
    else:
            selectionSort(myFidels, [startingRhos, totalCounts, myDensities, myfVals])
#
# checkSort([myFidels,startingRhos, totalCounts, myDensities, myfVals,perfectFidels,perfectDensities,perfectFVals,noCorrFidels,noCorrFVals,noCorrDensities],oldData)

#Print Website

if(any(settingsArray)):
    FORREPLACE = '<table>'
    FORREPLACE += '<tr><th>Settings used</th><th></th></tr>'
    FORREPLACE += '<tr><td><ul>'
    for x in range(len(settingsArray)):
        if settingsArray[x] :
            FORREPLACE += '<li>'
            FORREPLACE += settingsArrayName[x]
            FORREPLACE += '</li>'

    FORREPLACE += '</ul></td><td><img src="fidelGraph.png" style="width:80%;height:auto;"></td></tr></table>'
else:
    FORREPLACE = str('<img src="fidelGraph.png">')


for j in range(nStates):

    '''Real State'''
    FORREPLACE += '<table class="data">'
    FORREPLACE += '<tr><th colspan="2">Actual Densities</th><th>Counts</th>'
    # # Does not perform any tests.
    if(testPerfectTomo):
        FORREPLACE += '<th colspan="3">Perfect Densities</th>'

    FORREPLACE += '<th colspan="3">Caclulated Densities</th>'
    if (testCrossTalk):
        FORREPLACE += '<th colspan="2">CrossTalk</th>'
    FORREPLACE += '</tr>'
    FORREPLACE += '<tr>'
    FORREPLACE += '<td colspan="2">'
    FORREPLACE += printMatrix(startingRhos[j])
    FORREPLACE += '</td>'
    FORREPLACE += '<td>'
    FORREPLACE += '<div class="title">Counts: </div>' + str(int(AtotalCounts[j]))
    FORREPLACE += '</td>'

    # Does not perform any tests.
    if (testPerfectTomo):
        FORREPLACE += '<td>'
        FORREPLACE += printfValsFidelity(perfectFVals[j], perfectFidels[j])
        FORREPLACE += '</td>'
        FORREPLACE += '<td colspan="2">'
        FORREPLACE += printMatrix(perfectDensities[j])
        FORREPLACE += '</td>'

    '''My tomography'''
    FORREPLACE += '<td>'
    FORREPLACE += printfValsFidelity(myfVals[j],myFidels[j])
    FORREPLACE += '</td>'
    FORREPLACE += '<td colspan="2">'
    FORREPLACE += printMatrix(myDensities[j])
    FORREPLACE += '</td>'

    #Print crosstalk matrix
    if(testCrossTalk):
        FORREPLACE += '<td colspan="2">'
        FORREPLACE += printMatrix(myCTalks[j])
        FORREPLACE += '</td>'

    FORREPLACE += '</tr>'
FORREPLACE += '</table>'

fig = plt.figure()
fig.clf()

plt.plot(np.log10(AtotalCounts),myFidels,'.b')
if(testPerfectTomo and any(settingsArray)):
    plt.plot(np.log10(AtotalCounts), perfectFidels, '.g')
plt.title('Fidelities')
plt.xlabel("Log(Counts) base 10")
plt.ylabel("Fidelity")

plt.savefig("Results/Fidel_Graph_.png")

with open('Results/Template.html', 'r') as f:
    fff = '\n'.join(f.readlines())
    f.close()

fff = fff.replace('TOREPLACE', str(FORREPLACE))

with open('Results/DataVisualized'+str(numQubits)+'.html','w') as ff:
    ff.write(fff)
    ff.close()
print('Test Done. Results Printed')