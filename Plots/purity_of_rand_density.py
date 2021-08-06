import QuantumTomography as qLib
import numpy as np
import matplotlib.pyplot as plt


numSamples = 100000
numBins = 40



purities1 = np.zeros(numSamples)
purities2 = np.zeros(numSamples)
for i in range(numSamples):
    state1 = qLib.random_density_state(1)
    purities1[i] = qLib.purity(state1)
    state2 = qLib.random_density_state(2)
    purities2[i] = qLib.purity(state2)

print("Num Purities less then .5 : " + str(sum(purities1<.5)))
print("Num Purities less then .25 : " + str(sum(purities2<.25)))
[yVals1,binEdges1,foo] = plt.hist(purities1,bins=numBins,range=(.5,1))
[yVals2,binEdges2,foo] = plt.hist(purities2,bins=numBins,range=(.25,1))
yVals1 = yVals1/numSamples
yVals2 = yVals2/numSamples
xVals1 = np.zeros_like(yVals1)
xVals2 = np.zeros_like(yVals2)
for i in range(len(xVals1)):
    xVals1[i] = (binEdges1[i]+binEdges1[i+1])/2
for i in range(len(xVals2)):
    xVals2[i] = (binEdges2[i]+binEdges2[i+1])/2
plt.close()

fig, axes = plt.subplots(nrows=1, ncols=2, figsize=(9, 4))

# Plot 1 qubit
axes[0].scatter(xVals1, yVals1, s=20, facecolors='none', edgecolors='k')
axes[0].set_xlabel("Purity")
axes[0].set_ylabel("Frequency")
axes[0].set_title("Purity of random_density_state for 1 Qubit",pad=15)
axes[0].set_ylim((0,.08))

# Plot 2 qubit
axes[1].scatter(xVals2, yVals2, s=20, facecolors='none', edgecolors='k')
axes[1].set_xlabel("Purity")
axes[1].set_ylabel("Frequency")
axes[1].set_title("Purity of random_density_state for 2 Qubits",pad=15)
axes[1].set_ylim((0,.14))
plt.subplots_adjust(left=0.1,
                    bottom=0.15,
                    right=0.9,
                    top=0.85,
                    wspace=0.4,
                    hspace=0.4)

plt.savefig("images/purity_of_rand_density.png")