import unittest
from QuantumTomography.ProcessTomo import ProcessTomography
import numpy as np
import QuantumTomography as qlib
import matplotlib.pyplot as plt


class Test_ProcessTomo(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(Test_ProcessTomo, self).__init__(*args, **kwargs)

        self.states_6 = np.array([[1,0],[0,1],[0.7071,0.7071],[0.7071,-0.7071],[0.7071,0.7071j],[0.7071,-0.7071j]], dtype=complex)
        self.states_4 = np.array([[1,0],[0,1],[0.7071,0.7071],[0.7071,0.7071j]], dtype=complex)

        self.states_6_rhos = {}
        self.states_4_rhos = {}

        for i in range(6):
            self.states_6_rhos[i] = qlib.toDensity(self.states_6[i,:])
        for i in range(4):
            self.states_4_rhos[i] = qlib.toDensity(self.states_4[i,:])


    def test_pgates(self):
        num_measurements = 4
        num_steps = 100
        angles = []
        average_fidelities = []

        for k in range(num_steps):
            proportion = k / num_steps
            angle = proportion * (np.pi * 3)
            gate = np.array([[1, 0],[0, np.cos(angle) + 1j*np.sin(angle)]])
            output_actual = {}
            output_QPT = {}
            fidelities = []

            # getting the coincidences for 6 measurements
            if num_measurements == 6:
                coincidences = np.zeros((6,6))
                coincidences_predicted = np.zeros((6,6))
                for i in range(num_measurements):
                    output_actual[i] = qlib.toDensity(gate @ self.states_6[i,:])

                    for j in range(num_measurements):
                        coincidences[j,i] = 1000 * np.real(np.trace(self.states_6_rhos[j] @ output_actual[i]))
                chi = qlib.ProcessTomo.ProcessTomography(coincidences)

                for i in range(num_measurements):
                    output_QPT[i] = qlib.ProcessTomo.post_process_density(chi, self.states_6_rhos[i])
                    fidelities.append(qlib.fidelity(output_QPT[i], output_actual[i]))
                    for j in range(num_measurements):
                        coincidences_predicted[j,i] = 1000* np.real(np.trace(self.states_6_rhos[j] @ output_QPT[i]))


            # getting the coincidences for 4 measurements
            elif num_measurements == 4:
                coincidences = np.zeros((4,4))
                coincidences_predicted = np.zeros((4,4))
                for i in range(num_measurements):
                    output_actual[i] = qlib.toDensity(gate @ self.states_4[i,:])
                    for j in range(num_measurements):
                        coincidences[j,i] = 1000 * np.real(np.trace(self.states_4_rhos[j] @ output_actual[i]))
                chi = qlib.ProcessTomo.ProcessTomography(coincidences, self.states_4, self.states_4)
                for i in range(num_measurements):
                    output_QPT[i] = qlib.ProcessTomo.post_process_density(chi, self.states_4_rhos[i])
                    fidelities.append(qlib.fidelity(output_QPT[i], output_actual[i]))
                    for j in range(num_measurements):
                        coincidences_predicted[j, i] = 1000 * np.real(np.trace(self.states_4_rhos[j] @ output_QPT[i]))


            angles.append(angle)
            average_fidelities.append(np.average(fidelities))
        ax = plt.subplot(111)
        ax.plot(angles, average_fidelities)
        ax.ticklabel_format(useOffset=False)
        plt.show()


    def test_XYZH(self):
        average_fidelities_6 = []
        average_fidelities_4 = []
        gate_toplot = [0, 1, 2, 3]

        gates = {
            'X': np.array([[0, 1], [1, 0]], dtype=complex),
            'Y': np.array([[0, -1j], [1j, 0]], dtype=complex),
            'Z': np.array([[1, 0], [0, 1]], dtype=complex),
            'H': (1 / np.sqrt(2)) * np.array([[1, 1], [1, -1]], dtype=complex)
        }

        for k in range(len(gates.values())):
            gate = list(gates.values())[k]
            output_actual = {}
            output_QPT = {}
            fidelities_6 = []
            num_measurements = 6

            # getting the coincidences for 6 measurements
            if num_measurements == 6:
                coincidences = np.zeros((6,6))
                coincidences_predicted = np.zeros((6,6))
                for i in range(num_measurements):
                    output_actual[i] = qlib.toDensity(gate @ self.states_6[i,:])

                    for j in range(num_measurements):
                        coincidences[j,i] = 1000 * np.real(np.trace(self.states_6_rhos[j] @ output_actual[i]))
                chi = qlib.ProcessTomo.ProcessTomography(coincidences)

                for i in range(num_measurements):
                    output_QPT[i] = qlib.ProcessTomo.post_process_density(chi, self.states_6_rhos[i])
                    fidelities_6.append(qlib.fidelity(output_QPT[i], output_actual[i]))
                    for j in range(num_measurements):
                        coincidences_predicted[j,i] = 1000* np.real(np.trace(self.states_6_rhos[j] @ output_QPT[i]))

            average_fidelities_6.append(np.average(fidelities_6))


            num_measurements = 4
            fidelities_4 = []
            # getting the coincidences for 4 measurements
            if num_measurements == 4:
                coincidences = np.zeros((4,4))
                coincidences_predicted = np.zeros((4,4))
                for i in range(num_measurements):
                    output_actual[i] = qlib.toDensity(gate @ self.states_4[i,:])
                    for j in range(num_measurements):
                        coincidences[j,i] = 1000 * np.real(np.trace(self.states_4_rhos[j] @ output_actual[i]))
                chi = qlib.ProcessTomo.ProcessTomography(coincidences, self.states_4, self.states_4)
                for i in range(num_measurements):
                    output_QPT[i] = qlib.ProcessTomo.post_process_density(chi, self.states_4_rhos[i])
                    fidelities_4.append(qlib.fidelity(output_QPT[i], output_actual[i]))
                    for j in range(num_measurements):
                        coincidences_predicted[j, i] = 1000 * np.real(np.trace(self.states_4_rhos[j] @ output_QPT[i]))

            average_fidelities_4.append(np.average(fidelities_4))

        # ax6 = plt.subplot(121)
        # ax6.plot(gate_toplot, average_fidelities_6)
        # ax6.ticklabel_format(useOffset=False)
        #
        # ax4 = plt.subplot(122)
        # ax4.plot(gate_toplot, average_fidelities_4)
        # ax4.ticklabel_format(useOffset=False)
        # plt.show()






if __name__ == '__main__':
    unittest.main()