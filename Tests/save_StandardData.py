import numpy as np
import QuantumTomography as qk
import matplotlib.pyplot as plt
import scipy.io

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is used to run tomo against the standard_test_data, and save the Results"""

"WARNING! These tests run on the published library installed in your pip version, not the code in the local directory."

# Generate Fidelity Graph for the set standard testing data.
# Get a Chris's Data
mat = scipy.io.loadmat('Standard_Test_Data_1Q.mat')
mat2 = scipy.io.loadmat('Standard_Test_Data_2Q.mat')

# Set MLE Tomography
t1 = qk.Tomography()
t2 = qk.Tomography()
t1.setConfSetting("NQubits",1)
tomoInput1 = t1.getTomoInputTemplate()
tomoInput2 = t2.getTomoInputTemplate()

myFidels1_MLE = np.zeros_like(mat['fs'])
myFidels2_MLE = np.zeros_like(mat2['fs'])

# Set up the measurements
def getBasisMeas(numBits):
    basis = np.array([[1, 0], [0, 1], [(2 ** (-1 / 2)), (2 ** (-1 / 2))], [(2 ** (-1 / 2)), -(2 ** (-1 / 2))],
                      [(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j], [(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j]],
                     dtype=complex)
    m = basis.copy()
    for j in range(numBits-1):
        m = np.kron(m,basis)
    return m

Measurements_Density_1 = np.array([qk.toDensity(x) for x in getBasisMeas(1)])
Measurements_Density_2 = np.array([qk.toDensity(x) for x in getBasisMeas(2)])

for i in range(0,6):
    for j in range(0,100):
        # Put in values for 1qubit
        # Counts
        tomoInput1[:, 1 + 1] = mat['countsP'][i][j]
        for x in range(0,6):
            density = np.outer(tomoInput1[:, np.arange(1+2, 3*1+2)][x],tomoInput1[:, np.arange(1+2, 3*1+2)][x].conj())
            diff = density - mat["projOps"][x]
            if(any(abs(diff.flatten())>.0005)):
                print("STOP")
        # Put in values for 2qubit
        tomoInput2[:, 2 + 1] = mat2['countsP'][i][j]

        # ----------------------

        # MLE 1 qubit
        r1 = t1.state_tomography(tomoInput1)
        myFidels1_MLE[i][j] = qk.fidelity(mat["pReal"][i][j],r1[0])
        # MLE 2 qubit
        r2 = t2.state_tomography(tomoInput2)
        myFidels2_MLE[i][j] = qk.fidelity(mat2["pReal"][i][j], r2[0])

# Plot fidelities 1q
plt.title("Fidelities 1q")

CF = plt.scatter(np.log10(mat['totCounts']), mat['fs'], marker='x', color="b")
MF = plt.scatter(np.log10(mat['totCounts']), myFidels1_MLE, marker='x', color="r")

plt.legend((CF, MF),
           ('Chris', 'My MLE', 'My BAY'),
           scatterpoints=1,
           loc='lower right')
plt.ylabel("Fidelity")
plt.xlabel("Log(Total Counts)")

# plt.show()
plt.savefig('Results/Standard_test_data_q1.png')

# plotting fidelities 2q
plt.title("Fidelities 2q")
CF = plt.scatter(np.log10(mat2['totCounts']), mat2['fs'], marker='x', color="b")
MF = plt.scatter(np.log10(mat2['totCounts']), myFidels2_MLE, marker='x', color="r")

plt.legend((CF, MF),
           ('Chris', 'My MLE'),
           scatterpoints=1,
           loc='lower right')
plt.ylabel("Fidelity")
plt.xlabel("Log(Total Counts)")
# plt.show()
plt.savefig('Results/Standard_test_data_q2.png')

