import unittest
import QuantumTomography as qLib
import numpy as np
import numpy.testing as tests
from TestRun import runTests
import os

# Run tomographies
[[Tomo_Object_1, Fidelity_with_Original, Original_Purity, Total_Time]] = runTests(1,1)
[[Tomo_Object_2, Fidelity_with_Original, Original_Purity, Total_Time]] = runTests(2,1)
[[Tomo_Object_3, Fidelity_with_Original, Original_Purity, Total_Time]] = runTests(2,1,test2Det=True,testAccCorr=True,testCrossTalk=True,testDrift=True,testBell=True,errBounds=3,method="linear")


class Test_Functions(unittest.TestCase):

    def test_eval(self):
        filename = "Test_States/rand_eval.txt"
        for Tomo_Object in [Tomo_Object_1,Tomo_Object_2,Tomo_Object_3]:
            # Export
            Tomo_Object.exportToEval(filename)

            # Import
            Tomo_Object_copy = qLib.Tomography()
            Tomo_Object_copy.importEval(filename)

            # Make sure conf settings are the same
            for k in Tomo_Object_copy.conf.keys():
                if not (isinstance(Tomo_Object_copy.conf[k], np.ndarray)):
                    self.assertEqual(Tomo_Object.conf[k],Tomo_Object_copy.conf[k])

            # make sure inputs are the same
            tests.assert_array_equal(Tomo_Object.last_input,Tomo_Object_copy.last_input)
            tests.assert_array_equal(Tomo_Object.intensities,Tomo_Object_copy.intensities)

            # Run tomographies and make sure estimates are the same
            tests.assert_array_equal(Tomo_Object.last_rho, Tomo_Object_copy.last_rho)
            self.assertEqual(Tomo_Object.last_intensity, Tomo_Object_copy.last_intensity)
            self.assertEqual(Tomo_Object.last_fval, Tomo_Object_copy.last_fval)

    def test_conf_and_data(self):
        filename_c = "Test_States/rand_conf.txt"
        filename_d = "Test_States/rand_data.txt"
        for Tomo_Object in [Tomo_Object_1,Tomo_Object_2,Tomo_Object_3]:
            # Export
            Tomo_Object.exportToConf(filename_c)
            Tomo_Object.exportToData(filename_d)

            # Import
            Tomo_Object_copy = qLib.Tomography()
            Tomo_Object_copy.importConf(filename_c)
            Tomo_Object_copy.importData(filename_d)

            # Make sure conf settings are the same
            for k in Tomo_Object_copy.conf.keys():
                if not (isinstance(Tomo_Object_copy.conf[k], np.ndarray)):
                    self.assertEqual(Tomo_Object.conf[k],Tomo_Object_copy.conf[k])

            # make sure inputs are the same
            tests.assert_array_equal(Tomo_Object.last_input,Tomo_Object_copy.last_input)
            tests.assert_array_equal(Tomo_Object.intensities,Tomo_Object_copy.intensities)

            # Run tomographies and make sure estimates are the same
            tests.assert_array_equal(Tomo_Object.last_rho, Tomo_Object_copy.last_rho)
            self.assertEqual(Tomo_Object.last_intensity, Tomo_Object_copy.last_intensity)
            self.assertEqual(Tomo_Object.last_fval, Tomo_Object_copy.last_fval)

    def test_export_website(self):
        Fixed_Tomo_Object_1 = qLib.Tomography()
        Fixed_Tomo_Object_2 = qLib.Tomography()
        Fixed_Tomo_Object_3 = qLib.Tomography()
        Fixed_Tomo_Objs = [Fixed_Tomo_Object_1, Fixed_Tomo_Object_2, Fixed_Tomo_Object_3]

        for i in range(len(Fixed_Tomo_Objs)):
            Fixed_Tomo_Objs[i].importEval("Test_States/fixed_eval_"+str(i)+".txt")
            Fixed_Tomo_Objs[i].exportToConf_web("Test_States/Website_Files/conf_temp.txt")
            Fixed_Tomo_Objs[i].exportToData_web("Test_States/Website_Files/data_temp.txt")


            Fixed_Tomo_Objs[i].printLastOutput(bounds=10)
            print("----------------------")

            self.maxDiff = None
            self.assertEqual(open("Test_States/Website_Files/conf_"+str(i)+".txt").read(),
                             open("Test_States/Website_Files/conf_temp.txt").read())
            self.assertEqual(open("Test_States/Website_Files/data_" + str(i) + ".txt").read(),
                             open("Test_States/Website_Files/data_temp.txt").read())

    def test_printLastOutput(self):
        Fixed_Tomo_Object_1 = qLib.Tomography()
        Fixed_Tomo_Object_2 = qLib.Tomography()
        Fixed_Tomo_Object_3 = qLib.Tomography()
        Fixed_Tomo_Objs = [Fixed_Tomo_Object_1, Fixed_Tomo_Object_2, Fixed_Tomo_Object_3]

        for i in range(len(Fixed_Tomo_Objs)):
            Fixed_Tomo_Objs[i].importEval("Test_States/fixed_eval_" + str(i) + ".txt")
            Fixed_Tomo_Objs[i].printLastOutput()
            qLib.printLastOutput(Fixed_Tomo_Objs[i])


if __name__ == '__main__':
    unittest.main()
