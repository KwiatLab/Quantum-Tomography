import matplotlib.pyplot as plt
import numpy as np

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is used to graph the data created by the tests"""

maxNumberOfQubits = 1
timeVsCounts_xVals = []
timeVsCounts_yVals = []

for x in range(1,maxNumberOfQubits+1):

    '''### Fidelity Table for each qubit###'''

    # My data
    groupData = np.loadtxt('myLibrary'+str(x)+'.csv', delimiter=',')

    timeVsCounts_xVals = np.append(timeVsCounts_xVals,x*np.ones_like(groupData[2]))
    timeVsCounts_yVals = np.append(timeVsCounts_yVals,groupData[2])

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
plt.savefig("timeGraph_NQUBITS.png")
