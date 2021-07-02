import unittest
import QuantumTomography as qLib
import numpy as np
import numpy.testing as tests
from TestRun import runTests


class Test_Functions(unittest.TestCase):

    def __init__(self, *args, **kwargs):
        super(Test_Functions, self).__init__(*args, **kwargs)

        # initializing the inputs to be used in the following test tomographies
        self.tomoinp1 = np.array(
            [[1, 0, 500, 1, 0], [1, 0, 500, 0, 1], [1, 0, 500, 0.7071, 0.7071], [1, 0, 500, 0.7071, -0.7071],
             [1, 0, 500, 0.7071, 0.7071j], [1, 0, 500, 0.7071, -0.7071j]])
        self.tomoinp2 = np.array(
            [[1, 0, 0, 500, 1, 0, 1, 0], [1, 0, 0, 500, 1, 0, 0, 1], [1, 0, 0, 500, 1, 0, 0.7071, 0.7071],
             [1, 0, 0, 500, 1, 0, 0.7071, 0.7071j], [1, 0, 0, 500, 0, 1, 1, 0], [1, 0, 0, 500, 0, 1, 0, 1],
             [1, 0, 0, 500, 0, 1, 0.7071, 0.7071], [1, 0, 0, 500, 0, 1, 0.7071, 0.7071j],
             [1, 0, 0, 500, 0.7071, 0.7071, 1, 0], [1, 0, 0, 500, 0.7071, 0.7071, 0, 1],
             [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071], [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071j],
             [1, 0, 0, 500, 0.7071, 0.7071j, 1, 0], [1, 0, 0, 500, 0.7071, 0.7071j, 0, 1],
             [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071], [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071j]])

    def test_getCoincidences(self):
        # setting up tomography objects
        tomoObj1 = qLib.Tomography(1)
        tomoObj2 = qLib.Tomography(2)
        # setting up new inputs to be changed
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        # testing different random measurement values
        for i in range(500, 1500, 200):
            rand_coincidence1 = np.random.binomial(i, 0.5, 6)
            rand_coincidence2 = np.random.binomial(i, 0.5, 16)
            curinp1[:, 2] = rand_coincidence1
            curinp2[:, 3] = rand_coincidence2
            tomo = tomoObj1.StateTomography_Matrix(curinp1)
            tomo2 = tomoObj2.StateTomography_Matrix(curinp2)
            self.assertTrue((tomoObj1.getCoincidences() == rand_coincidence1).all())
            self.assertTrue((tomoObj2.getCoincidences() == rand_coincidence2).all())

    def test_getSingles(self):
        # setting up tomography objects
        tomoObj1 = qLib.Tomography(1)
        tomoObj2 = qLib.Tomography(2)
        # setting up new inputs to be changed
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        # test case with 0 everywhere
        tomo = tomoObj1.StateTomography_Matrix(curinp1)
        tomo2 = tomoObj2.StateTomography_Matrix(curinp2)
        self.assertTrue((tomoObj1.getSingles() == np.zeros(6)).all())
        self.assertTrue((tomoObj2.getSingles() == np.zeros((16, 2))).all())

        # test cases with different random singles values
        for i in range(1, 20, 4):
            singles_onebit = np.random.binomial(i, 0.5, 6)
            singles_twobit = np.random.binomial(i, 0.5, (16, 2))
            curinp1[:, 1] = singles_onebit
            curinp2[:, 1:3] = singles_twobit
            tomo = tomoObj1.StateTomography_Matrix(curinp1)
            tomo2 = tomoObj2.StateTomography_Matrix(curinp2)
            self.assertTrue((np.real(np.transpose(tomoObj1.getSingles()))[0] == singles_onebit).all())
            self.assertTrue((np.real(tomoObj2.getSingles()) == singles_twobit).all())

    def test_getTimes(self):
        # setting up tomography objects
        tomoObj1 = qLib.Tomography(1)
        tomoObj2 = qLib.Tomography(2)
        # make copies of inputs to change
        curinp1 = self.tomoinp1.copy()
        curinp2 = self.tomoinp2.copy()

        # test case with all ones for time
        tomo = tomoObj1.StateTomography_Matrix(curinp1)
        tomo2 = tomoObj2.StateTomography_Matrix(curinp2)
        self.assertTrue((tomoObj1.getTimes() == np.ones(6)).all())
        self.assertTrue((tomoObj2.getTimes() == np.ones(16)).all())

        # testing some different cases:
        for i in range(3):
            randTime1 = np.random.randint(1, 10, 6)
            randTime2 = np.random.randint(1, 10, 16)
            curinp1[:, 0] = randTime1
            curinp2[:, 0] = randTime2

            tomo = tomoObj1.StateTomography_Matrix(curinp1)
            tomo2 = tomoObj2.StateTomography_Matrix(curinp2)
            self.assertTrue((tomoObj1.getTimes() == randTime1).all())
            self.assertTrue((tomoObj2.getTimes() == randTime2).all())

    def test_buildTomoInput(self):
        tomoObj1 = qLib.Tomography(1)
        tomoObj2 = qLib.Tomography(2)

        [times1, singles1, coincidences1, measurements11, measurements12] = self.tomoinp1.transpose()

        measurements1 = np.array([measurements11, measurements12]).transpose()
        built_input1 = tomoObj1.buildTomoInput(measurements1, coincidences1, crosstalk=-1, efficiency=0, time=times1,
                                               singles=singles1, window=0, error=0)
        self.assertTrue((built_input1 == self.tomoinp1).all())

        [times2, singles21, singles22, coincidences2, measurements211, measurements212, measurements221,
         measurements222] = self.tomoinp2.transpose()
        singles2 = np.array([singles21, singles22]).transpose()

        measurements2 = np.array([measurements211, measurements212, measurements221, measurements222]).transpose()
        built_input2 = tomoObj2.buildTomoInput(measurements2, coincidences2, crosstalk=-1, efficiency=0, time=times2,
                                               singles=singles2, window=0, error=0)

        self.assertTrue((built_input2 == self.tomoinp2).all())


    # This tests phaserToComplex,complexToPhaser,removeGlobalPhase and stateToString
    def test_stateToString(self):
        phase1 = np.array([0.786236, 3.366])
        complex1 = qLib.phaserToComplex(phase1)
        phase1_2 = qLib.complexToPhaser(complex1)
        complex1_2 = qLib.phaserToComplex(phase1_2)
        tests.assert_array_almost_equal(phase1,phase1_2)
        tests.assert_almost_equal(complex1,complex1_2)

        complex2 = .5 + 5j
        phase2 = qLib.complexToPhaser(complex2)
        complex2_2 = qLib.phaserToComplex(phase2)
        phase2_2 = qLib.complexToPhaser(complex2_2)
        tests.assert_array_almost_equal(phase2,phase2_2)
        tests.assert_almost_equal(complex2,complex2_2)

        state = np.array([-.5 + .5j, -.5 - .5j])

        self.assertEqual('0.707|H> + i0.707|V>',qLib.stateToString(state))

    # Error bounds is not specified, then getProperties(5) is run after
    # then do getProperties(3) is run and should not increase the number of monte carlo states.
    def test_getProperties(self):
        for i in [1, 2]:
            tomo = runTests(numQubits=i, nStates=1)[0][0]
            self.assertEqual(len(tomo.mont_carlo_states), 1)
            tomo.getProperties(5)
            self.assertEqual(len(tomo.mont_carlo_states), 6)
            tomo.getProperties(3)
            self.assertEqual(len(tomo.mont_carlo_states), 6)

    def test_printLastOutput(self):
        for i in [1,2]:
            tomo = runTests(numQubits=i, nStates=1)[0][0]
            self.assertTrue(True)

    # This tests several functions involving going from density to tvals
    # functions include t_to_density(), t_matrix(), density2tm(), and density2t()
    def test_tvals(self):
        for n in range(1,5):
            og_state = qLib.random_density_state(n)
            curr_state = og_state
            for i in range(100):
                tvals = qLib.density2t(curr_state)
                next_state = qLib.t_to_density(tvals)
                fid = qLib.fidelity(next_state,og_state)
                self.assertGreater(fid,.95)
                curr_state = next_state


if __name__ == '__main__':
    unittest.main()
