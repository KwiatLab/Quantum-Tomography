import numpy as np
import matplotlib.pyplot as plt

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

"""This script is used to save data to the results folder. Create an instance of the class and use addData.
This class will automatically save the data after a certain amount of new entires are given"""

class csvSaver():
    def __init__(self,numQubits):
        self.n = numQubits
        self.previousData = self.loadData()
        self.count = 0

        # determines how many states it should hold before saving it to the csv file.
        saveNumber = ["Error",500,100,5]
        self.nStates = 1
        try:
            self.nStates = saveNumber[self.n]
        except :
            self.nStates = 1
        self.newData = np.zeros((3,self.nStates))

    def loadData(self):
        try:
            return np.loadtxt('Results/myLibrary'+str(self.n)+'.csv', delimiter=',')
        except:
            return np.zeros((3,0))
    def saveData(self):
        finalData = np.concatenate((self.previousData,self.newData[:, ~np.all(self.newData == 0, axis=0)]), axis=1)
        np.savetxt('Results/myLibrary'+str(self.n)+'.csv',finalData, delimiter=',')
        self.previousData = finalData
        self.newData = np.zeros((3, self.nStates))
    def addData(self,count,fidel,time):
        if(fidel!=-1):
            self.newData[:,self.count] = [count,fidel,time]
            self.count +=1
            if(self.count>= self.nStates):
                self.saveData()
                self.count = 0
