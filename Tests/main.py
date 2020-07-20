"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

from __future__ import print_function
__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE : http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

import argparse
import os
import matplotlib.pyplot as plt
import warnings
warnings.filterwarnings('ignore')
from TestRun import TestRun
from PIL import Image

def file_path(string):
    if os.path.isfile(string) or string == "../../" :
        return string
    else:
        raise FileNotFoundError(string)

def dir_path(string):
    if os.path.isdir(string) or string == "../../" :
        return string
    else:
        raise FileNotFoundError(string)

def main():
    # create argument parser object
    parser = argparse.ArgumentParser(description="Quantum Tomography")

    parser.add_argument("-f", "--full", action='store_true',default=False,
                        help="Do a full test with all possible setting combinations. Default Tests each settings seperate")
    parser.add_argument("-w", "--web", action='store_true',default=False,
                        help="Including this will make the tests go on the web")

    # parse the arguments from standard input
    args = parser.parse_args()

    doFullTest = args.full
    useWebsite = args.web
    numFails = 0

    class bcolors:
        OKBLUE = '\033[94m'
        OKGREEN = '\033[92m'
        FAIL = '\033[91m'

    print(f"{bcolors.OKBLUE}STARTING TESTS")

    if doFullTest:
        #  Test([nBits, nStat, det2, cross, bounds, bell, acc, drift, web])=
        for nbits in range(1,4):
            for det2 in range(0, 2):
                for cross in range(0, 2):
                    for bounds in [0,3]:
                        for bell in range(0, 2):
                            for acc in range(0, 2):
                                for drift in range(0, 2):
                                    numFails+= TestRun([nbits, bounds, acc,det2, cross, bell, drift, useWebsite]).run()
    else:
         #TestRun([nBits, bounds, acc, det2, cross, bell, drift, web])
        numFails += TestRun([1, 0, 0, 0, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([1, 3, 0, 0, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([1, 0, 1, 0, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([1, 0, 0, 1, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([1, 0, 0, 0, 1, 0, 0, useWebsite]).run()
        numFails += TestRun([1, 0, 0, 0, 0, 1, 0, useWebsite]).run()
        numFails += TestRun([1, 0, 0, 0, 0, 0, 1, useWebsite]).run()
        numFails += TestRun([1, 0, 1, 0, 0, 0, 1, useWebsite]).run()

        numFails += TestRun([2, 0, 0, 0, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([2, 3, 0, 0, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([2, 0, 1, 0, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([2, 0, 0, 1, 0, 0, 0, useWebsite]).run()
        numFails += TestRun([2, 0, 0, 0, 1, 0, 0, useWebsite]).run()
        numFails += TestRun([2, 0, 0, 0, 0, 1, 0, useWebsite]).run()
        numFails += TestRun([2, 0, 0, 0, 0, 0, 1, useWebsite]).run()

    print(f"{bcolors.OKBLUE}ALL TESTS COMPLETED")
    fig = plt.figure()
    fig.clf()
    plt.plot(TestRun.allCounts, TestRun.allFidels, '.b')
    plt.title('Fidelities')
    plt.xlabel("Log(Counts) base 10")
    plt.ylabel("Fidelity")
    plt.savefig('Results/bigGraph.png')

    if (numFails == 0):
        print(f"{bcolors.OKGREEN}ALL TESTS PASSED! [" + str(TestRun.counter-numFails)+"/" +str(TestRun.counter)+"]")
    else:
        print(f"{bcolors.FAIL}[" + str(numFails)+"/" +str(TestRun.counter)+"] TESTS FAILED")

    im = Image.open(r"C:\Users\scoot\Desktop\Kwiat\Python_Testing\Testing\Results\bigGraph.png")
    im.show()
if __name__ == "__main__":
    main()