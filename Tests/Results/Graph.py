import matplotlib.pyplot as plt
import numpy as np
from QuantumTomography import floatToString


maxNumberOfQubits = 1
timeVsCounts_xVals = []
timeVsCounts_yVals = []
# For descriptive stats at bottom of image
p = np.array([['nBits',"nStates",'Mean',"std"],['nQ',"nStates",'Mean',"std"],['nQ',"nStates",'Mean',"std"],['nQ',"nStates",'Mean',"std"]],dtype="O")

for x in range(1,maxNumberOfQubits+1):

    '''### Fidelity Table for each qubit###'''

    # My data
    groupData = np.loadtxt('myLibrary'+str(x)+'.csv', delimiter=',')

    timeVsCounts_xVals = np.append(timeVsCounts_xVals,x*np.ones_like(groupData[2]))
    timeVsCounts_yVals = np.append(timeVsCounts_yVals,groupData[2])
    p[x, 0] = x
    p[x, 1] = len(groupData[2])
    p[x, 2] = np.mean(groupData[2])
    p[x, 3] = np.std(groupData[2])
    fig = plt.figure()
    fig.clf()
    plt.plot(np.log10(groupData[0]),groupData[1],'.b',label="Our Group's Tomography")
    #
    # # Other Source
    # plt.plot(np.log10(AtotalCounts), perfectFidels, '.g')
    plt.legend(loc='lower left')
    plt.title('Fidelities')
    plt.xlabel("Log(Counts) base 10")
    plt.ylabel("Fidelity")
    plt.savefig("fidelGraph"+str(x)+".png")

    '''### TIME vs Counts Table ###'''
    fig = plt.figure()
    fig.clf()
    plt.plot(np.log10(groupData[0]), groupData[2], '.b',label="Our Group's Tomography")
    plt.legend(loc='upper left')
    plt.title('Time')
    plt.xlabel("Log(Counts) base 10")
    plt.ylabel("Time (sec)")
    plt.savefig("timeGraph" + str(x) + ".png")

'''### Time vs Nqubits Graph ###'''


fig = plt.figure()
fig.clf()
plt.plot(timeVsCounts_xVals, timeVsCounts_yVals, '.b',label="Our Group's Tomography")
plt.legend(loc='upper left')
plt.title('Time')
plt.xlabel("nQubits")
plt.ylabel("Time (sec)")
fig.subplots_adjust(bottom=0.3)
# Adds descriptive stats at bottom of image

mx = 0
descriptiveStats=""
for i in range(p.shape[0]):
    for j in range(p.shape[1]):
        if(not isinstance(p[i, j],str)):
            p[i, j] = floatToString(p[i, j]).replace(" ", "") + "  "
        if(len(p[i, j])>mx):
            mx = len(p[i, j])

for i in range(p.shape[0]):
    for j in range(p.shape[1]):
        descriptiveStats += p[i, j] + " "*(mx-len(p[i, j]))
    descriptiveStats += "\n"


plt.figtext(0.2, 0.02, descriptiveStats, wrap=True, horizontalalignment='left', fontsize=12, fontname='monospace')

plt.savefig("timeGraph_NQUBITS.png")
