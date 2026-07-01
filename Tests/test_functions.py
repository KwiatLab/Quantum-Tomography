import numpy as np
import numpy.testing as tests
import pytest
from TestRun import runTests

import QuantumTomography as qLib


@pytest.fixture
def tomoObj1():
    return qLib.Tomography(1)


@pytest.fixture
def tomoObj2():
    return qLib.Tomography(2)


@pytest.fixture
def tomoinp1():
    return np.array(
        [
            [1, 0, 500, 1, 0],
            [1, 0, 500, 0, 1],
            [1, 0, 500, 0.7071, 0.7071],
            [1, 0, 500, 0.7071, -0.7071],
            [1, 0, 500, 0.7071, 0.7071j],
            [1, 0, 500, 0.7071, -0.7071j],
        ]
    )


@pytest.fixture
def tomoinp2():
    return np.array(
        [
            [1, 0, 0, 500, 1, 0, 1, 0],
            [1, 0, 0, 500, 1, 0, 0, 1],
            [1, 0, 0, 500, 1, 0, 0.7071, 0.7071],
            [1, 0, 0, 500, 1, 0, 0.7071, 0.7071j],
            [1, 0, 0, 500, 0, 1, 1, 0],
            [1, 0, 0, 500, 0, 1, 0, 1],
            [1, 0, 0, 500, 0, 1, 0.7071, 0.7071],
            [1, 0, 0, 500, 0, 1, 0.7071, 0.7071j],
            [1, 0, 0, 500, 0.7071, 0.7071, 1, 0],
            [1, 0, 0, 500, 0.7071, 0.7071, 0, 1],
            [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071],
            [1, 0, 0, 500, 0.7071, 0.7071, 0.7071, 0.7071j],
            [1, 0, 0, 500, 0.7071, 0.7071j, 1, 0],
            [1, 0, 0, 500, 0.7071, 0.7071j, 0, 1],
            [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071],
            [1, 0, 0, 500, 0.7071, 0.7071j, 0.7071, 0.7071j],
        ]
    )


def test_getCoincidences(tomoinp1, tomoinp2, tomoObj1, tomoObj2):

    # testing different random measurement values
    for i in range(500, 1500, 200):
        rand_coincidence1 = np.random.binomial(i, 0.5, 6)
        rand_coincidence2 = np.random.binomial(i, 0.5, 16)
        tomoinp1[:, 2] = rand_coincidence1
        tomoinp2[:, 3] = rand_coincidence2
        tomo = tomoObj1.StateTomography_Matrix(tomoinp1)
        tomo2 = tomoObj2.StateTomography_Matrix(tomoinp2)
        assert np.array_equal(tomoObj1.getCoincidences(), rand_coincidence1)
        assert np.array_equal(tomoObj2.getCoincidences(), rand_coincidence2)


def test_getSingles(tomoinp1, tomoinp2, tomoObj1, tomoObj2):
    # test case with 0 everywhere
    tomo = tomoObj1.StateTomography_Matrix(tomoinp1)
    tomo2 = tomoObj2.StateTomography_Matrix(tomoinp2)
    assert np.array_equal(tomoObj1.getSingles(), np.zeros((6, 1)))
    assert np.array_equal(tomoObj2.getSingles(), np.zeros((16, 2)))

    # test cases with different random singles values
    for i in range(1, 20, 4):
        singles_onebit = np.random.binomial(i, 0.5, 6)
        singles_twobit = np.random.binomial(i, 0.5, (16, 2))
        tomoinp1[:, 1] = singles_onebit
        tomoinp2[:, 1:3] = singles_twobit
        tomo = tomoObj1.StateTomography_Matrix(tomoinp1)
        tomo2 = tomoObj2.StateTomography_Matrix(tomoinp2)
        assert np.array_equal(np.real(np.transpose(tomoObj1.getSingles()))[0], singles_onebit)
        assert np.array_equal(np.real(tomoObj2.getSingles()), singles_twobit)


def test_getTimes(tomoinp1, tomoinp2, tomoObj1, tomoObj2):
    # test case with all ones for time
    tomo = tomoObj1.StateTomography_Matrix(tomoinp1)
    tomo2 = tomoObj2.StateTomography_Matrix(tomoinp2)
    assert np.array_equal(tomoObj1.getTimes(), np.ones(6))
    assert np.array_equal(tomoObj2.getTimes(), np.ones(16))

    # testing some different cases:
    for i in range(3):
        randTime1 = np.random.randint(1, 10, 6)
        randTime2 = np.random.randint(1, 10, 16)
        tomoinp1[:, 0] = randTime1
        tomoinp2[:, 0] = randTime2

        tomo = tomoObj1.StateTomography_Matrix(tomoinp1)
        tomo2 = tomoObj2.StateTomography_Matrix(tomoinp2)
        assert np.array_equal(tomoObj1.getTimes(), randTime1)
        assert np.array_equal(tomoObj2.getTimes(), randTime2)


def test_buildTomoInput(tomoinp1, tomoinp2, tomoObj1, tomoObj2):

    [times1, singles1, coincidences1, measurements11, measurements12] = tomoinp1.transpose()

    measurements1 = np.array([measurements11, measurements12]).transpose()
    built_input1 = tomoObj1.buildTomoInput(
        measurements1, coincidences1, crosstalk=-1, efficiency=0, time=times1, singles=singles1, window=0, error=0
    )
    assert np.array_equal(built_input1, tomoinp1)

    [
        times2,
        singles21,
        singles22,
        coincidences2,
        measurements211,
        measurements212,
        measurements221,
        measurements222,
    ] = tomoinp2.transpose()
    singles2 = np.array([singles21, singles22]).transpose()

    measurements2 = np.array([measurements211, measurements212, measurements221, measurements222]).transpose()
    built_input2 = tomoObj2.buildTomoInput(
        measurements2, coincidences2, crosstalk=-1, efficiency=0, time=times2, singles=singles2, window=0, error=0
    )

    assert np.array_equal(built_input2, tomoinp2)


# This tests phaserToComplex,complexToPhaser,removeGlobalPhase and stateToString
def test_stateToString(tomoinp1, tomoinp2, tomoObj1, tomoObj2):
    phase1 = np.array([0.786236, 3.366])
    complex1 = qLib.phaserToComplex(phase1)
    phase1_2 = qLib.complexToPhaser(complex1)
    complex1_2 = qLib.phaserToComplex(phase1_2)
    tests.assert_array_almost_equal(phase1, phase1_2)
    tests.assert_almost_equal(complex1, complex1_2)

    complex2 = 0.5 + 5j
    phase2 = qLib.complexToPhaser(complex2)
    complex2_2 = qLib.phaserToComplex(phase2)
    phase2_2 = qLib.complexToPhaser(complex2_2)
    tests.assert_array_almost_equal(phase2, phase2_2)
    tests.assert_almost_equal(complex2, complex2_2)

    state = np.array([-0.5 + 0.5j, -0.5 - 0.5j])

    assert qLib.stateToString(state) == "0.707|H> + i0.707|V>"


# Error bounds is not specified, then getProperties(5) is run after
# then do getProperties(3) is run and should not increase the number of monte carlo states.
def test_getProperties(tomoinp1, tomoinp2, tomoObj1, tomoObj2):
    for i in [1, 2]:
        tomo = runTests(numQubits=i, nStates=1)[0][0]
        assert len(tomo.mont_carlo_states) == 1
        tomo.getProperties(5)
        assert len(tomo.mont_carlo_states) == 6
        tomo.getProperties(3)
        assert len(tomo.mont_carlo_states) == 6


# This tests several functions involving going from density to tvals
# functions include t_to_density(), t_matrix(), density2tm(), and density2t()
def test_tvals(tomoinp1, tomoinp2, tomoObj1, tomoObj2):
    for n in range(1, 5):
        og_state = qLib.random_density_state(n)
        curr_state = og_state
        for i in range(100):
            tvals = qLib.density2t(curr_state)
            next_state = qLib.t_to_density(tvals)
            fid = qLib.fidelity(next_state, og_state)
            assert fid > 0.95
            curr_state = next_state
