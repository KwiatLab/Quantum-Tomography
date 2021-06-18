from __future__ import print_function
from .TomoFunctions import *
from .TomoClassHelpers import *
from .Utilities import ConfDict,getValidFileName
import numpy as np
from scipy.optimize import leastsq,minimize
from scipy.linalg import cholesky
import warnings

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""


"""
    class Tomography()
    Desc: This is the main tomography object that the library is built around. The goal is to only have one tomography object and edit the
    configuration settings as you go.
    """
class Tomography():

    # # # # # # # # # # # #
    '''Public Variables'''
    # # # # # # # # # # # #

    # conf['NQubits']: >= 1, it will take a much longer time for more qubits.
    # conf['NDetectors']: 1 or 2
    # conf['Crosstalk']: [[C0->0, C0->1], [C1->0, C1->1]]
    # conf['Bellstate']: 'no' or 'yes'
    # conf['DoDriftCorrection'] = 'no' or 'yes'
    # conf['DoAccidentalCorrection'] = 'no' or 'yes'
    # conf['DoErrorEstimation']: >= 0
    # conf['Window']: 0 or array like, dimension = 1
    # conf['Efficiency']: 0 or array like, dimension = 1
    # conf['Beta']: 0 to 0.5, depending on purity of state and total number of measurements.

    # tomo_input: array like, dimension = 2.
    # input data of the last tomography run.

    # For n detectors:
    # tomo_input[:, 0]: times
    # tomo_input[:, np.arange(1, n_qubit+1)]: singles
    # tomo_input[:, n_qubit+1]: coincidences
    # tomo_input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
    #
    # For 2n detectors:
    # tomo_input[:, 0]: times
    # tomo_input[:, np.arange(1, 2*n_qubit+1)]: singles
    # tomo_input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
    # tomo_input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements

    # last_rho: The predicted density matrix of the last tomography run.
    # last_intensity: The predicted intensity of the state for the last tomography run.
    # last_fval: The final value of the internal optimization function for the last tomography run.
    # last_input: The last tomo_input matrix used
    # mont_carlo_states: The generated monte carlo states.
    # intensity: Relative pump power (arb. units) during measurement for the last tomography run.
    # err_functions: These are the properties that you want to be calculated when getProperties is called.

    # # # # # # # # # #
    '''Constructors'''
    # # # # # # # # # #

    """
    Default Constructor
    Desc: This initializes a default tomography object
    """
    def __init__(self,nQ = 2):
        self.conf = ConfDict([('NQubits',nQ), ('NDetectors', 1),
                              ('Crosstalk', np.array([[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]])),
                              ('Bellstate', 0),('DoDriftCorrection', 0), ('DoAccidentalCorrection', 0),
                              ('DoErrorEstimation', 0),('Window', 0),('Efficiency', 0),('RhoStart', []),
                              ('Beta', 0)])
        self.err_functions = ['concurrence', 'tangle', 'entropy', 'linear_entropy', 'negativity', 'purity']

    """
    setConfSetting(setting, val)
    Desc: Sets a specific self.conf setting. This is no longer needed but will be kept in the class to
    support old code.

    Parameters
    ----------
    setting : string
        The setting you want to update.
        Possible values are ['NQubits', 'NDetectors', 'Crosstalk', 'Bellstate', 'DoDriftCorrection', 'DoAccidentalCorrection', 'DoErrorEstimation', 'Window', 'Efficiency', 'RhoStart', 'Beta']
    val: ndarray, int, or string
            The new value you want to the setting to be.
    """
    def setConfSetting(self, setting, val):
        raise DeprecationWarning('As of v1.0.3.7 setConfSetting() is no longer needed to set a conf setting. '
                                 'Settings can now be set directly like a normal dictionary.')
        if (isinstance(val, str)):
            if (valC.lower() == "yes" or val.lower() == "true"):
                valC = 1
            elif (valC.lower() == "no" or val.lower() == "false"):
                valC = 0
        self.conf[setting] = valC

    """
    importConf(conftxt)
    Desc: Import a text file containing the configuration settings.

    Parameters
    ----------
    conftxt : string
        path to configuration file
    """
    def importConf(self, conftxt):
        conf = self.conf
        exec(compile(open(conftxt, "rb").read(), conftxt, 'exec'))

    """
    importData(datatxt)
    Desc: Import a text file containing the tomography data and run tomography.

    Parameters
    ----------
    datatxt : string
       path to data file
    Returns
    -------
    rhog : ndarray with shape = (2^numQubits, 2^numQubits)
        The predicted density matrix.
    intensity : float
        The predicted overall intensity used to normalize the state.
    fvalp : float
        Final value of the internal optimization function. Values greater than the number
        of measurements indicate poor agreement with a quantum state.
    """
    def importData(self, datatxt):
        exec(compile(open(datatxt, "rb").read(), datatxt, 'exec'))
        return self.state_tomography(locals().get('tomo_input'), locals().get('intensity'))

    """
    importEval(evaltxt)
    Desc: Import a eval file containing the tomography data and the configuration setting, then run tomography.

    Parameters
    ----------
    evaltxt : string
        path to eval file
    Returns
    -------
    rhog : ndarray with shape = (2^numQubits, 2^numQubits)
        The predicted density matrix.
    intensity : float
        The predicted overall intensity used to normalize the state.
    fvalp : float
        Final value of the internal optimization function. Values greater than the number
        of measurements indicate poor agreement with a quantum state.
    """
    def importEval(self, evaltxt):
        conf = self.conf
        exec(compile(open(evaltxt, "rb").read(), evaltxt, 'exec'))
        return self.state_tomography(locals().get('tomo_input'), locals().get('intensity'))

    # # # # # # # # # # # # # #
    '''Tomography Functions'''
    # # # # # # # # # # # # # #

    """
    state_tomography(tomo_input, intensities)
    Desc: Main function that runs tomography.
    Parameters
    ----------
    tomo_input : ndarray
        The input data for the current tomography. This is what tomo_input will be set to. Example can be seen at top of page.
        See getTomoInputTemplate() to get a template for this input.
    intensities : 1darray with length = number of measurements
        Relative pump power (arb. units) during measurement; used for drift correction. Default will be an array of ones
    Returns
    -------
    rhog : ndarray with shape = (2^numQubits, 2^numQubits)
        The predicted density matrix.
    intensity : The predicted overall intensity used to normalize the state.
        The predicted overall intensity used to normalize the state.
    fvalp : float
        Final value of the internal optimization function. Values greater than the number
        of measurements indicate poor agreement with a quantum state.
    """
    def state_tomography(self, tomo_input, intensities = -1,method="MLE",saveState=True):

        # define a uniform intenstiy if not stated
        if (isinstance(intensities, int)):
            self.intensities = np.ones(tomo_input.shape[0])
        elif(len(intensities.shape) == 1 and intensities.shape[0] == tomo_input.shape[0]):
            self.intensities = intensities
        else:
            raise ValueError("Invalid intensities array")

        # filter the data
        self.last_input = tomo_input
        [coincidences, measurements_densities, measurements_pures, accidentals] = self.filter_data(tomo_input)
        self.conf['Method'] = method

        # get the starting state from linear_tomography if not defined
        starting_matrix = self.conf['RhoStart']

        try:
            # todo: go over linear tomography. Clean it up. Figure out why it fails on the very rare occasion.
            if not starting_matrix:
                starting_matrix = self.linear_tomography(coincidences, measurements_pures)[0]
            # Currently linear tomography gets the phase wrong. So a temporary fix is to just transpose it.
            starting_matrix = starting_matrix.transpose()
            starting_matrix = make_positive(starting_matrix)
            starting_matrix = starting_matrix / np.trace(starting_matrix)
        except:
            raise RuntimeError('Failed to run linear Tomography')

        if(method == "MLE"):
            # perform MLE tomography
            [rhog, intensity, fvalp] = self.tomography_MLE(starting_matrix, coincidences, measurements_densities, accidentals)
        else:
            # perform Bayesian tomography
            [rhog, intensity, fvalp] = self.tomography_BME(starting_matrix, coincidences, measurements_densities,accidentals)

        if saveState:
            # save the results
            self.last_rho = rhog.copy()
            self.last_intensity = intensity
            self.last_fval = fvalp
            # Reset the monte carlo states after doing a tomography
            self.mont_carlo_states = list([[rhog, intensity, fvalp]])

        return [rhog, intensity, fvalp]

    """
        state_tomo(measurements, counts)
        todo:comment
        """
    def state_tomo(self,measurements,counts,crosstalk=-1,efficiency=0,time=-1,singles=-1,window=0,error=0, intensities=-1,method="MLE"):
        tomo_input = self.buildTomoInput(measurements, counts, crosstalk, efficiency, time, singles,window,error)
        return self.state_tomography(tomo_input, intensities, method=method)



    """
    tomography_MLE(starting_matrix, coincidences, m, acc)
    Desc: Calculates the most likely state given the data.

    Parameters
    ----------
    starting_matrix : ndarray with shape = (2^numQubits, 2^numQubits)
        The starting predicted state found with linear tomography.
    coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The counts of the tomography.
    m : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
        The measurements of the tomography in density matrix form.
    acc : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.
    Returns
    -------
    rhog : ndarray with shape = (2^numQubits, 2^numQubits)
        The predicted density matrix.
    intensity : float
        The predicted overall intensity used to normalize the state.
    fvalp : float
        Final value of the internal optimization function. Values greater than the number
        of measurements indicate poor agreement with a quantum state.
    """
    def tomography_MLE(self, starting_matrix, coincidences, m, accidentals):
        init_intensity = np.mean(np.multiply(coincidences, 1/self.intensities)) * (starting_matrix.shape[0])
        starting_tvals = density2t(starting_matrix)
        n_t = len(starting_tvals)
        # np.multiply(coincidences, self.intensities)
        # t_to_density(starting_tvals)
        starting_tvals = starting_tvals + 0.0001
        starting_tvals = starting_tvals * np.sqrt(init_intensity)

        coincidences = np.real(coincidences)
        coincidences = coincidences.flatten()

        bet = self.conf['Beta']
        # n_data = np.shape(coincidences)[0]


        prediction = np.zeros(m.shape[2]) + 0j

        if bet == 0:
            final_tvals = leastsq(self.maxlike_fitness, np.real(starting_tvals), args = (coincidences, accidentals, m, prediction))[0]
            fvalp = np.sum(self.maxlike_fitness(final_tvals, coincidences, accidentals, m, prediction) ** 2)
        else:
            final_tvals = \
            leastsq(self.maxlike_fitness_hedged, np.real(starting_tvals), args = (coincidences, accidentals, m, prediction, bet))[0]
            fvalp = np.sum(self.maxlike_fitness_hedged(final_tvals, coincidences, accidentals, m, prediction, bet) ** 2)

        final_matrix = t_to_density(final_tvals)
        intensity = np.trace(final_matrix)
        final_matrix = final_matrix / np.trace(final_matrix)

        intensity = np.float64(np.real(intensity))

        return [final_matrix, intensity, fvalp]

    """
    maxlike_fitness(t, coincidences, accidentals, m, prediction)
    Desc: Calculates the diffrence between the current predicted state data and the actual data.

    Parameters
    ----------
    t : ndarray
        T values of the current predicted state.
    coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The counts of the tomography.
    accidentals : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.
    m : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
        The measurements of the tomography in density matrix form.
    prediction : ndarray
        Predicted counts from the predicted state.

    Returns
    -------
    val : float
        value of the optimization function.
    """
    def maxlike_fitness(self, t, coincidences, accidentals, m, prediction):
        rhog = t_to_density(t)
        for j in range(len(prediction)):
            prediction[j] = np.float64(np.real(self.intensities[j] * np.real(np.trace(np.dot(m[:, :, j], rhog))) + accidentals[j]))
            prediction[j] = np.max([prediction[j], 0.01])
            # Avoid dividing by zero for pure states
            if (prediction[j] == 0):
                prediction[j] == 1

        val = (prediction - coincidences) / np.sqrt(prediction)

        val = np.float64(np.real(val))

        return val


    """
    maxlike_fitness_hedged(t, coincidences, accidentals, m, prediction, bet)
    Desc: Calculates the diffrence between the current predicted state data and the actual data using hedged maximum likelihood.

    Parameters
    ----------
    t : ndarray
        T values of the current predicted state.
    coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The counts of the tomography.
    accidentals : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.
    m : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
        The measurements of the tomography in density matrix form .
    prediction : ndarray
        Predicted counts from the predicted state.
    bet : float 0 to .5
        The beta value used.

    Returns
    -------
    val : float
        value of the optimization function.
    """
    def maxlike_fitness_hedged(self, t, coincidences, accidentals, m, prediction, bet):

        rhog = t_to_density(t)

        for j in range(len(prediction)):
            prediction[j] = self.intensities[j] * np.real(np.trace(np.dot(m[:, :, j], rhog))) + accidentals[j]
            prediction[j] = np.max([prediction[j], 0.01])

        hedge = np.repeat(np.real((bet * np.log(np.linalg.det(np.mat(rhog)))) / len(prediction)), len(prediction))
        val = np.sqrt(np.real((((prediction - coincidences) ** 2) / (2 * prediction)) - hedge) + 1000)

        val = np.float64(np.real(val))

        return val

    def tomography_BME(self,starting_matrix, counts, measurements, accidentals):

        # algorithm parameters
        preallocationSize = 100000  # preallocate memory for speed
        changeThreshold = 10 ** -9  # stopping condition threshold
        numberOfPriorSamples = 5000  # number of samples to sample from the inital prior
        updateFrequency = 1000  # controls how often the estimate of the posterior is updated, and how often stability is checked
        max_iterations = 150  # The maximum number of times we update the posterior

        # initialize and preallocatememory
        d = starting_matrix.shape[0]
        sectionLikelihood = 0
        cumLikelihood = 0
        i = 0
        samplesToSC = 0
        iterationCounter = 0
        eps = 2.2204e-16

        # Lists that are evaluated periodically
        stabilityHistory = np.zeros(max_iterations+1)
        # stoppingConditionTestLocations = np.zeros_like(stabilityHistory)

        # List of Sample States
        sampleStates_rho = np.zeros((d, d, preallocationSize), dtype=complex)
        sampleStates_loglike = np.zeros(preallocationSize)
        sampleStates_tvals = np.zeros((preallocationSize, d ** 2))

        # Current best guess for posterior parameters
        mean_tvals = np.zeros(d ** 2)

        # estimate the state intensity for the guess state
        norm = sum(counts / (len(counts) * d))
        minimization_output = minimize(self.log_likelyhood, norm,
                              args=(density2t(starting_matrix), counts, measurements, accidentals))
        norm = minimization_output.x[0]
        base_loglike = minimization_output.fun

        # SAMPLE FROM PRIOR
        # Samples are drawn from the ginibre distribution
        stillUsingPrior = 1
        while stillUsingPrior:
            # expand memory in chunks if preallocation size is exceeded
            if (i >= len(sampleStates_loglike)):
                sampleStates_rho = np.concatenate((sampleStates_rho, np.zeros((d, d, i))), axis=2)
                sampleStates_loglike = np.concatenate((sampleStates_loglike, np.zeros(i)))
                sampleStates_tvals = np.concatenate((sampleStates_tvals, np.zeros((i, d ** 2))), axis=0)

            # Sample from the ginibre distribution
            random_rho = random_density_state(self.conf['NQubits'])
            random_tvals = density2t(random_rho)
            random_loglike = self.log_likelyhood(norm, random_rho, counts, measurements, accidentals)

            # Store sample state, tval, and its loglikelyhood
            sampleStates_rho[:, :, i] = random_rho
            sampleStates_tvals[i, :] = random_tvals
            sampleStates_loglike[i] = random_loglike - base_loglike

            # If we have met the min required amount of prior samples try to switching to the posterior estimate
            if i >= numberOfPriorSamples:
                # if the range of likelihoods becomes toolarge, rescale
                if min(sampleStates_loglike[:i+1])<0:
                    base_loglike = min(sampleStates_loglike[:i+1])
                    sampleStates_loglike[:i+1] = sampleStates_loglike[:i+1] - base_loglike

                    # find the most optimal intensity for this best state
                    minimization_output = minimize(self.log_likelyhood, norm,
                                          args=(random_rho, counts, measurements, accidentals))
                    norm = minimization_output.x[0]

                # Calculate the covariance matrix for the posterior
                sampleStates_like = normalizeExponentLikelihood(sampleStates_loglike[:i+1])
                contributingParams = sampleStates_like > 0
                firstMax_like = max(sampleStates_like)
                secondMax_like = max(sampleStates_like[sampleStates_like < firstMax_like])
                [mean_tvals, cov_tvals] = weightedcov(sampleStates_tvals[:i + 1][contributingParams, :],
                                                           sampleStates_like[contributingParams])

                # switch to adaptive distributions only if the covariance matrix is valid
                # cov must be positive semi - definite
                try:
                    if not all(np.linalg.eigvals(cov_tvals)>=0) or (secondMax_like / firstMax_like) < .01:
                        # Keep using the prior
                        stillUsingPrior = 1
                        numberOfPriorSamples = numberOfPriorSamples + updateFrequency
                    else:
                        # Stop using the prior
                        stillUsingPrior = 0
                except:
                    # Keep using the prior
                    stillUsingPrior = 1
                    numberOfPriorSamples = numberOfPriorSamples + updateFrequency

            # track the total likelihood and the section likelihood
            # cumLikelihood=cumLikelihood+np.exp(-1*sampleStates_loglike[i])
            # sectionLikelihood = sectionLikelihood + np.exp(-1*sampleStates_loglike[i])

            # Increase iteration number
            i = i + 1

        # Calculate the running mean and likelihood with the states sampled from the prior
        unNormedMean_rho = np.dot(sampleStates_rho[:, :, :i], np.exp(-1 * sampleStates_loglike[:i]))
        tr = np.trace(unNormedMean_rho)
        if tr == 0:
            mean_rho = np.zeros(d)
        else:
            mean_rho = unNormedMean_rho / tr
        prevMean_rho = mean_rho

        # SAMPLE FROM POSTERIOR
        # Samples by drawing tVals from a multivariate normal distro
        stillUsingPosterior = 1
        while stillUsingPosterior:
            # expand memory in chunks if preallocation size is exceeded
            if (i >= len(sampleStates_loglike)):
                sampleStates_rho = np.concatenate((sampleStates_rho, np.zeros((d, d, i))), axis=2)
                sampleStates_loglike = np.concatenate((sampleStates_loglike, np.zeros(i)))
                sampleStates_tvals = np.concatenate((sampleStates_tvals, np.zeros((i, d ** 2))), axis=0)

            # Draws by sampling the tVals from a multivariate normal distro
            random_tvals = rand.multivariate_normal(mean_tvals, cov_tvals)
            random_rho = t_to_density(random_tvals)
            random_loglike = self.log_likelyhood(norm, random_rho, counts, measurements, accidentals)

            # Store sample state, tval, and its likelyhood
            sampleStates_loglike[i] = random_loglike - base_loglike
            sampleStates_rho[:, :, i] = random_rho
            sampleStates_tvals[i, :] = random_tvals

            # update the running mean and likelihood with the newly sampled state
            unNormedMean_rho = unNormedMean_rho + sampleStates_rho[:, :, i] * np.exp(-1 * sampleStates_loglike[i])
            cumLikelihood = cumLikelihood + np.exp(-1 * sampleStates_loglike[i])
            sectionLikelihood = sectionLikelihood + np.exp(-1 * sampleStates_loglike[i])

            # Periodically perform the following tasks
            if i % updateFrequency == 0:

                # Update the posterior parameters
                # -------------------------------

                # if the range of likelihoods becomes toolarge, rescale
                if min(sampleStates_loglike[:i + 1]) < 0:
                    base_loglike = min(sampleStates_loglike[:i + 1])
                    sampleStates_loglike[:i + 1] = sampleStates_loglike[:i + 1] - base_loglike
                    unNormedMean_rho = np.dot(sampleStates_rho[:, :, :i], np.exp(-1 * sampleStates_loglike[:i]))

                    # find the most optimal intensity for this best state
                    minimization_output = minimize(self.log_likelyhood, norm,
                                                   args=(random_rho, counts, measurements, accidentals))
                    norm = minimization_output.x[0]

                sampleStates_like = normalizeExponentLikelihood(sampleStates_loglike[0:i + 1])
                contributingParams = sampleStates_like > 0
                firstMax_like = max(sampleStates_like)
                secondMax_like = max(sampleStates_like[sampleStates_like < firstMax_like])
                [newMean_tvals, newCov_tvals] = weightedcov(sampleStates_tvals[:i + 1][contributingParams, :],
                                                                 sampleStates_like[contributingParams])

                # Update the posterior estimate only if the covariance matrix is valid
                # cov must be positive semidefinite
                try:
                    if all(np.linalg.eigvals(cov_tvals) >= 0) and (secondMax_like / firstMax_like) > .001:
                        cov_tvals = newCov_tvals
                        mean_tvals = newMean_tvals
                except:
                    pass

                # Check if the stopping condition is met
                # --------------------------------------

                tr = np.trace(unNormedMean_rho)
                if tr == 0:
                    mean_rho = np.zeros(d)
                    stabilityHistory[iterationCounter] = float('inf')
                    # stoppingConditionTestLocations[iterationCounter] = i
                else:
                    mean_rho = unNormedMean_rho / tr

                    # calculate the normalization factor based on the fraction
                    # of total likelihood sampled in this iteration
                    normFactor = np.floor(i / updateFrequency) * sectionLikelihood / cumLikelihood
                    infidelity = (1 - np.real(fidelity(mean_rho, prevMean_rho)))
                    # ensure rounding never leads the infidelity to be negative
                    if infidelity < 0:
                        infidelity = eps

                    stabilityHistory[iterationCounter] = infidelity / normFactor
                    # stoppingConditionTestLocations[iterationCounter] = i

                    # Check if this is the first time we've calculated the posterior estimate
                    if (iterationCounter > 0) and (stabilityHistory[iterationCounter] < changeThreshold):
                        samplesToSC = i
                        stillUsingPosterior = 0

                    # IMPOSE STOPPING CONDITION IF TOO MANY ITERATIONS REACHED
                    if (iterationCounter > max_iterations):
                        warnings.warn("MAX ITER REACHED - BME. Stability: " + str(stabilityHistory[iterationCounter]))
                        samplesToSC = i
                        stillUsingPosterior = 0

                    # if we haven't reached stability, set the state to compare to next
                    # iteration, and reset the section likelihood
                    iterationCounter = iterationCounter + 1
                    prevMean_rho = mean_rho
                    sectionLikelihood = 0

            # Increase iteration number
            i = i + 1
        return [mean_rho, norm, samplesToSC]

    # todo: comment block
    def likelyhood(self,givenState,coincidences,measurments,accidentals,intensities):

        Averages = np.zeros_like(coincidences)

        for j in range(len(Averages)):
            Averages[j] =  intensities[j]*np.trace(np.matmul(measurments[:, :, j],givenState)) + accidentals[j]
            Averages[j] = np.max([Averages[j], 0.0000001])

        val = (Averages - coincidences)**2 / (2*Averages)
        val = np.float64(np.real(val))
        prob = np.exp(-1*np.sum(val),dtype=np.longdouble)
        return prob

    def log_likelyhood(self, intensity, givenState, coincidences, measurements, accidentals):

        if (len(givenState.shape) == 1):
            givenState = t_to_density(givenState)

        # # todo : currently the intensity is sometimes being passed in as a 1x1 array
        # try:
        #     intensity = intensity[0]
        # except:
        #     pass

        Averages = np.zeros_like(coincidences, dtype=np.float)
        for j in range(len(Averages)):
            Averages[j] = intensity * self.intensities[j] * np.real(np.trace(np.matmul(measurements[:, :, j], givenState))) + \
                          accidentals[j]

            # Avoid dividing by zero for pure states
            if (Averages[j] == 0):
                Averages[j] == 1

        val = (Averages - coincidences) ** 2 / (2 * Averages)
        return sum(val)

    """
    linear_tomography(coincidences, measurements)
    Desc: Uses linear techniques to find a starting state for maximum likelihood estimation.
    TODO: idk what this is based off of. Does it even do the correct thing?
    Parameters
    ----------
    coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The counts of the tomography.
    measurements : ndarray with shape = (number of measurements, 2*numQubits)
        The measurements of the tomography in pure state form .

    Returns
    -------
    rhog : ndarray with shape = (2^numQubits, 2^numQubits)
        The starting predicted state.
    intensity : float
        The predicted overall intensity used to normalize the state.
    """
    def linear_tomography(self, coincidences, measurements, m_set = ()):
        if m_set == ():
            m_set = independent_set(measurements)
        if np.isscalar(m_set):
            n = len(coincidences)
            linear_measurements = measurements
            linear_data = coincidences

        else:
            n = np.int(np.sum(m_set))
            linear_measurements = measurements[(np.rot90(m_set == 1.0)[0])]
            linear_data = coincidences[(np.rot90(m_set == 1.0)[0])]

        linear_rhog = np.zeros([measurements.shape[1], measurements.shape[1]])

        b = b_matrix(linear_measurements)
        b_inv = np.linalg.inv(b)

        m = np.zeros([measurements.shape[1], measurements.shape[1], n]) + 0j
        for j in range(n):
            m[:, :, j] = m_matrix(j, linear_measurements, b_inv)
            linear_rhog = linear_rhog + linear_data[j] * m[:, :, j]

        intensity = np.trace(linear_rhog)
        rhog = linear_rhog / intensity

        return [rhog, intensity]

    """
    filter_data(tomo_input, intensities)
    Desc: Filters the data into separate arrays.
    TODO: edit intensities out of comment

    Parameters
    ----------
    tomo_input : ndarray
        The input data for the current tomography. This is what self.last_input will be set to. Example can be seen at top of page.
        See getTomoInputTemplate() to get a template for this input.
    intensities : 1darray with length = number of measurements
        Relative pump power (arb. units) during measurement; used for drift correction.

    Returns
    -------
    data : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The counts/coincidences of the tomography.
    measurements_densities : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
        The measurements of the tomography in density matrix form.
    measurements_pures : ndarray with shape = (number of measurements, 2*numQubits)
        The measurements of the tomography in pure state form.
    acc : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.
    """
    def filter_data(self, tomo_input):
        # getting variables

        nbits = self.conf['NQubits']
        ndet = self.conf['NDetectors']

        # time values
        t = tomo_input[:, 0]

        # singles values
        n_singles = self.getNumSingles()
        # sings = tomo_input[:, np.arange(n_singles) + 1]
        sings = self.getSingles()

        # coincidences
        n_coinc = self.getNumCoinc()
        coinc = self.getCoincidences()

        # settings
        settings = tomo_input[:, np.arange(n_singles + n_coinc + 1, len(tomo_input[0]))]

        # Set window is not already defined
        if (isinstance(self.conf['Window'], int)):
            self.conf['Window'] = np.ones(n_coinc)
        eff = self.conf['Window'][0:n_coinc]

        # Accidental Correction
        acc = np.zeros_like(coinc)
        if(self.conf['DoAccidentalCorrection'] == 1):
            window = self.conf['Window']
            if np.isscalar(window):
                window = window*np.ones(n_coinc)
            scalerIndex = np.concatenate((np.ones(nbits - 2), [2, 2]))
            additiveIndex = np.array([0, 1])
            for j in range(2, nbits):
                additiveIndex = np.concatenate(([2*j], additiveIndex))
            if (len(coinc.shape) == 1):
                acc = acc[:, np.newaxis]
            for j in range(n_coinc):
                index = bin(j).split("b")[1]
                index = "0" * (nbits - len(index)) + index
                index = [int(char) for char in index]
                index = index*scalerIndex + additiveIndex
                index = np.array(index, dtype = int)
                acc[:, j] = np.prod(np.real(sings[:, tuple(index)]), axis = 1) * (window[j] * 1e-9 / np.real(t)) ** (nbits - 1)
            if (acc.shape != coinc.shape):
                acc = acc[:, 0]

        # crosstalk
        ctalk = np.array(self.conf['Crosstalk'])[0:2 ** nbits, 0:2 ** nbits]

        # This chunk of code handles edge cases of the input of the crosstalk matrix.
        # The important crosstalk matrix is the big_crosstalk
        # Todo : this could use some clean up
        crosstalk = ctalk
        if np.ndim(ctalk) >= 3:
            for j in range(ctalk.shape[2]):
                crosstalk[j] = ctalk[:, :, j]
        big_crosstalk = crosstalk[:]
        big_crosstalk = big_crosstalk * np.outer(eff, np.ones(n_coinc))

        # Get measurements
        # Todo: this is very messy but works...? At some point it can be cleaned up, so that toDensity is used
        # and the shapes are similar. The axis of the measuremennts_densities is kinda wonky. Not similar to measurements_pures
        measurements_densities = np.zeros([2 ** nbits, 2 ** nbits, np.prod(coinc.shape)],dtype=complex)
        measurements_pures = np.zeros([np.prod(coinc.shape), 2 ** nbits],dtype=complex)
        for j in range(coinc.shape[0]):
            m_twiddle = np.zeros([2 ** nbits, 2 ** nbits, 2 ** nbits]) + 0j
            u = 1
            for k in range(nbits):
                alpha = settings[j][2 * k]
                beta = settings[j][2 * k + 1]
                psi_k = np.array([alpha, beta])
                psip_k = np.array([np.conj(beta), np.conj(-alpha)])
                u_k = np.outer((np.array([1, 0])), psi_k) + np.outer((np.array([0, 1])), psip_k)
                # Changed from tensor_product to np.kron
                u = np.kron(u, u_k)
            if (ndet == 1):
                for k in range(0, 2 ** nbits):
                    m_twiddle[k, :, :] = np.outer(u[:, k].conj().transpose(), u[:, k])

                measurements_pures[j * n_coinc, :] = u[:, 0].conj().transpose()
                for k in range(1):
                    for l in range(2 ** nbits):
                        measurements_densities[:, :, j + k] = measurements_densities[:, :, j + k] + m_twiddle[l, :, :] * big_crosstalk[k, l]
                coincidences = coinc

            else:
                for k in range(2 ** nbits):
                    m_twiddle[k, :, :] = np.outer(u[:, k].conj().transpose(), u[:, k])
                    measurements_pures[j * n_coinc + k, :] = u[:, k].conj().transpose()
                for k in range(2 ** nbits):
                    for l in range(2 ** nbits):
                        measurements_densities[:, :, j * (2 ** nbits) + k] = measurements_densities[:, :, j * (2 ** nbits) + k] + m_twiddle[l, :, :] * \
                                                        big_crosstalk[k, l]

                coincidences = coinc.reshape((np.prod(coinc.shape), 1))
                acc = acc.reshape((np.prod(acc.shape), 1))


        # If 2det/qubit then expand intensity array
        if(ndet==2):
            self.intensities = np.kron(self.intensities,np.ones(2**nbits))

        return [coincidences, measurements_densities, measurements_pures, acc]

    """
        buildTomoInput(tomo_input, intensities)
        TODO: Test and write comment for this"""
    def buildTomoInput(self, measurements, counts, crosstalk, efficiency,time,singles,window, error):
        ################
        # measurements #
        ################
        # Get the number of qubits based on the measurement matrix
        if (len(measurements.shape)!=2 or measurements.shape[1]%2 != 0):
            raise ValueError("Invalid measurements matrix")
        else:
            self.conf['NQubits'] = int(measurements.shape[1]/2)
        # TODO: make sure measurements span the entire space

        ##########
        # counts #
        ##########
        # Check if counts has right dimensions
        if (counts.shape[0] != measurements.shape[0]):
            raise ValueError("Number of Counts does not match the number of measurements")
        # Determine if 2det/qubit from counts matrix
        try:
            if(counts.shape[1] == 1):
                self.conf['NDetectors'] = 1
            elif(counts.shape[1] == self.conf['NQubits']*2):
                self.conf['NDetectors'] = 2
            else:
                raise ValueError("The second axis of counts does not have the right dimension. Should be 1 or 2*NQubits for 2det")
        except:
            self.conf['NDetectors'] = 1

        ##############
        # efficiency #
        ##############
        if((not isinstance(efficiency, int)) and (len(efficiency.shape) !=1 or efficiency.shape[0] != self.conf['NQubits']*2)):
            raise ValueError("Invalid efficiency array. Length should be NQubits*2")
        else:
            self.conf['Efficiency'] = efficiency

        ##########
        # window #
        ##########
        if ((not isinstance(window, int)) and (len(window.shape) != 1 or window.shape[0] != self.getNumCoinc())):
            raise ValueError("Invalid window array. Length should be NQubits*NDetectors")
        else:
            self.conf['Window'] = window

        ################
        # time/singles #
        ################
        if ((not isinstance(time, int)) or (not isinstance(singles, int)) or (not isinstance(window, int))):
            self.conf['DoAccidentalCorrection'] = 1
            # Check if time has right length
            if (isinstance(time, int)):
                time = np.ones(measurements.shape[0])
            elif not (len(time.shape) == 1 and time.shape[0] == measurements.shape[0]):
                ValueError("Invalid time array")
            # Check if singles has right dimensions
            if (isinstance(singles, int)):
                singles = np.zeros((measurements.shape[0],2*self.conf["NQubits"]))
            elif not (len(singles.shape) == 2 and singles.shape[0] == measurements.shape[0] and singles.shape[1] == 2*self.conf["NQubits"]):
                raise ValueError("Invalid singles matrix")
            # Check if window has right length
            if (isinstance(window, int)):
                self.conf['Window'] = window
            elif (len(window.shape) == 1 and window.shape[0] == self.getNumCoinc()):
                self.conf['Window'] = window
            else:
                raise ValueError("Invalid window array")
        else:
            time = np.ones(measurements.shape[0])
            singles = np.zeros((measurements.shape[0],self.conf["NQubits"]))
            self.conf['Window'] = window


        #############
        # crosstalk #
        #############
        if (isinstance(crosstalk, int)):
            self.conf['Crosstalk'] = np.identity(2 ** self.conf["NQubits"])
        elif (len(crosstalk.shape) == 2 and crosstalk.shape[0] == crosstalk.shape[1] and crosstalk.shape[0] == 2 ** self.conf["NQubits"] ):
            self.conf['Crosstalk'] =  crosstalk
        else:
            raise ValueError("Invalid crosstalk matrix")

        #########
        # error #
        #########
        self.conf['DoErrorEstimation'] = error

        ##############
        # tomo_input #
        ##############
        # Here we build the tomo_input matrix and then return this to be fed to state_tomography
        tomo_input = 0
        n_qubit = self.conf['NQubits']
        if (self.conf['NDetectors'] == 1):
            tomo_input = np.zeros((measurements.shape[0], 3 * n_qubit + 2), dtype=complex)
            # times
            tomo_input[:, 0] = time
            # singles
            tomo_input[:, 1: n_qubit + 1] = singles
            # counts
            tomo_input[:, n_qubit + 1] = counts
            # measurements
            tomo_input[:, n_qubit + 2: 3 * n_qubit + 2] = measurements
        else:
            tomo_input = np.zeros((measurements.shape[0], 2 ** n_qubit + 4 * n_qubit + 1), dtype=complex)
            # times
            tomo_input[:, 0] = times
            # singles
            tomo_input[:, 1: 2 * n_qubit + 1] = singles
            # counts
            tomo_input[:, 2 * n_qubit + 1: 2 ** n_qubit + 2 * n_qubit + 1] = counts
            # measurements
            tomo_input[:, 2 ** n_qubit + 2 * n_qubit + 1: 2 ** n_qubit + 4 * n_qubit + 1] = measurements


        return tomo_input

    # # # # # # # # # #
    '''Get Functions'''
    # # # # # # # # # #


    """
    getNumCoinc()
    Desc: Returns the number of coincidences per measurement for the current configurations.
    """
    def getNumCoinc(self):
        if (self.conf['NDetectors'] == 2):
            return 2 ** self.conf['NQubits']
        else:
            return 1

    """
    getNumSingles()
    Desc: Returns the number of singles per measurement for the current configurations.
    """
    def getNumSingles(self):
        if (self.conf['NDetectors'] == 2):
            return 2 * self.conf['NQubits']
        else:
           return self.conf['NQubits']

    """
    getCoincidences()
    Desc: Returns an array of counts for all the measurments.
    """
    def getCoincidences(self):
        if (self.conf['NDetectors'] == 2):
            return np.real(self.last_input[:, np.arange(2*self.conf['NQubits']+1, 2**self.conf['NQubits']+2*self.conf['NQubits']+1)])
        else:
            return np.real(self.last_input[:, self.conf['NQubits']+1])

    """
    getSingles()
    Desc: Returns an array of singles for all the measurments.
    """
    def getSingles(self):
        if (self.conf['NDetectors'] == 2):
            return self.last_input[:, np.arange(1, 2*self.conf['NQubits']+1)]
        else:
            return self.last_input[:, np.arange(1, self.conf['NQubits']+1)]

    """
    getTimes()
    Desc: Returns an array of times for all the measurments.
    """
    def getTimes(self):
        return self.last_input[:, 0]

    """
        getMeasurements()
        Desc: Returns an array of measurements in pure state form for all the measurments.
        """
    def getMeasurements(self):
        if (self.conf['NDetectors'] == 2):
            return self.last_input[:, np.arange(2**self.conf['NQubits']+2*self.conf['NQubits']+1, 2**self.conf['NQubits']+4*self.conf['NQubits']+1)]
        else:
            return self.last_input[:, np.arange(self.conf['NQubits']+2, 3*self.conf['NQubits']+2)]

    """
    getNumBits()
    Desc: returns the number of qubits for the current configurations.
    """
    def getNumBits(self):
        return self.conf['NQubits']

    """
    getNumDetPerQubit()
    Desc: returns the number of detectors per qubit for the current configurations.
    """
    def getNumDetPerQubit(self):
        return self.conf['NDetectors']

    """
    getNumOfDetectorsTotal()
    Desc: Returns the total number of detectors for the current configurations.
    """
    def getNumOfDetectorsTotal(self):
        return self.getNumDetPerQubit()*self.getNumBits()

    """
    numBits(numBits)
    Desc: Returns an array of standard measurments in pure state form for the given number of qubits.
    TODO: edit comment to add numDet
    Parameters
    ----------
    numBits : int
        number of qubits you want for each measurement. Default will use the number of qubits in the current configurations.
    """
    def getBasisMeas(self, numBits = -1,numDet = -1):
        # check if numBits is an int
        if not (isinstance(numBits, int)):
            raise ValueError("numBits must be an integer")
        if(numBits < 1):
            numBits = self.getNumBits()

        # check if numDet is not 1 or 2
        if not (numDet in [-1, 1, 2]):
            raise ValueError("numDet must be an 1 or 2")
        if (numDet < 1):
            numDet = self.getNumDetPerQubit()

        # Define start meas basis
        if(numDet == 1):
            basis = np.array([[1, 0],
                              [0, 1],
                              [(2 ** (-1 / 2)), (2 ** (-1 / 2))],
                              [(2 ** (-1 / 2)), -(2 ** (-1 / 2))],
                              [(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j],
                              [(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j]],dtype=complex)
        else:
            basis = np.array([[1, 0],
                              [(2 ** (-1 / 2)), (2 ** (-1 / 2))],
                              [(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j]], dtype=complex)
        numBasis = basis.shape[0]
        m = np.zeros((numBasis ** numBits, 2 * numBits), dtype = complex)
        for i in range(m.shape[0]):
            for j in range(0, m.shape[1], 2):
                bitNumber = np.floor((m.shape[1] - j - 1) / 2)
                index = int(((i) / numBasis ** (bitNumber)) % numBasis)

                m[i, j] = basis[index][0]
                m[i, j + 1] = basis[index][1]
        return m

    """
    getTomoInputTemplate()
    Desc: returns a standard template for tomo_input for the given number of qubits.
    TODO: edit this comment for numDet
    Parameters
    ----------
    numBits : int
        number of qubits you want for each measurement. Default will use the number of qubits in the current configurations.

    Returns
    ----------
    Tomoinput : ndarray
        The input data for the current tomography. This is what tomo_input will be set to. Example can be seen at top of page.
        See getTomoInputTemplate() to get a template for this input.
    """
    def getTomoInputTemplate(self, numBits = -1,numDet = -1):
        # check if numBits is an int
        if not (isinstance(numBits, int)):
            raise ValueError("numBits must be an integer")
        if (numBits < 1):
            numBits = self.getNumBits()

        # check if numDet is not 1 or 2
        if not(numDet in [-1,1,2]):
            raise ValueError("numDet must be an 1 or 2")
        if (numDet < 1):
            numDet = self.getNumDetPerQubit()

        measurements = self.getBasisMeas(numBits,numDet)

        if(numDet == 1):
            # For n detectors:
            Tomoinput = np.zeros((6**numBits, 3*numBits+2), dtype = complex)

            # input[:, n_qubit+1]: coincidences
            # coincidences left zeros

            # input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
            Tomoinput[:, np.arange(numBits + 2, 3 * numBits + 2)] = measurements

        else:
            # For 2n detectors:
            Tomoinput = np.zeros((3**numBits, 2**numBits+4*numBits+1), dtype = complex)

            # input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
            # coincidences left zeros

            # input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements
            Tomoinput[:, np.arange(2**numBits+2*numBits+1, 2**numBits+4*numBits+1)] = measurements

        # tomo_input[:, 0]: times
        Tomoinput[:, 0] = np.ones_like(Tomoinput[:, 0])

        # tomo_input[:, np.arange(1, 2 * n_qubit + 1)]: singles
        # singles left zero
        return Tomoinput

    # # # # # # # # # # #
    '''ERROR FUNCTIONS'''
    # # # # # # # # # # #

    """
    getProperties(rho, bounds)
    Desc: Returns all the properties of the given density matrix.
          Using bounds will not change the conf settings. The calculated properties are determined by self.err_functions.

    Parameters
    ----------
    bounds : boolean
        Set this to true if you want error bounds on your estimated property values. Default is False.
        These are determined with monte carlo simulation and the states are saved under self.mont_carlo_states

    Returns
    -------
    vals : ndarray with shape = (length of self.err_functions, 2)
        The first col is the name of the property.
        The second col is the value of the property.
        The third col is the error bound on the property.
    """
    def getProperties(self,bounds = -1):
        # if bounds not set use the conf settings
        if bounds == -1:
            bounds = self.conf['DoErrorEstimation']

        # generate states if needed
        if bounds > len(self.mont_carlo_states)-1:
            self.tomography_states_generator(bounds-len(self.mont_carlo_states)+1)
        states = np.array(self.mont_carlo_states)[:bounds+1,:]
        if states.shape[0] == 1:
            vals1 = np.array([["intensity",np.mean(states[:,1]),"NA"],
                              ["fval", np.mean(states[:,2]),"NA"]])
        else:
            vals1 = np.array([["intensity", np.mean(states[:, 1]), np.std(states[:, 1],ddof=1)],
                              ["fval", np.mean(states[:, 2]), np.std(states[:, 2],ddof=1)]])
        vals2 = getProperties_helper_bounds(self.err_functions, states[:,0])

        vals = np.concatenate((vals1, vals2))
        return vals



    """
    tomography_states_generator(n)
    Desc: Uses monte carlo simulation to create random states similar to the estimated state.
          The states are also saved under self.mont_carlo_states

    Parameters
    ----------
    n : int
        Number of approximate states you want to create.
        If no value is given it will use the conf settings.
    Returns
    -------
    rhop : ndarray with shape = (n, 2^numQubits, 2^numQubits)
        The approximate density matrices.
    intenp : ndarray with length = n
        The intensity of each approximate density matrices.
    fvalp : ndarray with length = n
        The fval of the regression associated with each approximate density matrices.
    """
    def tomography_states_generator(self, n):
        # Save the last data so we can restore it later
        last_input = self.last_input.copy()

        ndet = self.conf['NDetectors']
        acc = self.conf['DoAccidentalCorrection']
        time = self.getTimes()
        meas = self.getMeasurements()
        singles = self.getSingles()
        counts = self.getCoincidences()
        for j in range(n):
            if ndet == 1:
                test_counts = np.zeros_like(counts)
                for k in range(len(counts)):
                        test_counts[k] = np.random.poisson(counts[k])
                if acc:
                    test_singles = np.zeros_like(singles)
                    for k in range(len(singles)):
                        test_singles[k] = np.random.poisson(singles[k])
                else:
                    test_singles = singles
                test_data = np.concatenate((np.array([time]).T, test_singles,np.array([test_counts]).T, meas), axis = 1)
                [rhop, intenp, fvalp] = self.state_tomography(test_data, self.intensities,saveState=False)
            elif ndet == 2:
                test_counts = np.zeros_like(counts)
                for k in range(len(counts)):
                    test_counts[k] = np.random.poisson(counts[k])
                if acc:
                    test_singles = np.zeros_like(singles)
                    for k in range(len(singles)):
                        test_singles[k] = np.random.poisson(singles[k])
                else:
                    test_singles = singles
                test_data = np.concatenate((np.array([time]).T, test_singles, np.array([test_counts]).T, meas), axis=1)
                [rhop, intenp, fvalp] = self.state_tomography(test_data, self.intensities,saveState=False)
            self.mont_carlo_states.append([rhop, intenp, fvalp])
        # Restore the last tomo_input matrix to the original one
        self.last_input = last_input
        return self.mont_carlo_states

    """
        getBellSettings(rho, bounds = -1)
        Desc: Returns the optimal measurment settings for the CHSH bell inequality.
              Using bounds will not change the conf settings.

        DISCLAIMER : In Progress, have not checked.

        Parameters
        ----------
        rho : ndarray with shape = (2^numQubits, 2^numQubits)
            The density matrix you would like to know the optimal bell measurment settings of.
        bounds : boolean
            Set this to true if you want error bounds on your estimated measurment settings. Default will use the conf settings.
            These are determined with monte carlo simulation and the states are saved under self.mont_carlo_states

        Returns
        -------
        vals : ndarray with shape = (length of self.err_functions, 2 or 3)
            The first col is the associated side
            The second col is the detector settings on the associated side.
            The third col is the error bound on the detector settings.
        """
    def getBellSettings(self, rho, partsize_init = 9, partsize = 5, t = 3, bounds = -1):
        # if bounds not set use the conf settings
        if (type(bounds) == int):
            bounds = (self.conf['DoErrorEstimation'] > 1) or self.mont_carlo_states != 0
        if (bounds):
            # check if given state is of the current data
            if (fidelity(self.last_rho, rho) < .95):
                raise ValueError("Input data needed. You can only calculate bounds on the state used in the tomography.")
            err_n = self.conf['DoErrorEstimation']
            # increase err_n if needed
            if (err_n <= 1):
                err_n = 2
            # generate states if needed
            if (self.mont_carlo_states == 0):
                states = self.tomography_states_generator(err_n)
            else:
                states = self.mont_carlo_states
                err_n = len(states[1])
            vals = getBellSettings_helper_bounds(states[0], rho, partsize_init, partsize, t, err_n)
        else:
            vals = getBellSettings_helper(rho, partsize_init, partsize, t)
        return vals

    """
        printLastOutput(tomo, bounds)
        Desc: Prints the properties of the last tomography to the console. Properties are defined in tomography conf settings.
              Using bounds will not change the conf settings. The calculated properties are determined by self.err_functions.

        Parameters
        ----------
        bounds : boolean
            The amount of monte carlo states to use in the error estimation. Default will use whatever is set in the conf settings.
            The states are saved under self.mont_carlo_states
        """
    def printLastOutput(self, bounds = -1):
        p = np.array(self.last_rho.copy(), dtype = "O")
        print("State: ")
        mx = 0
        for i in range(p.shape[0]):
            for j in range(p.shape[1]):
                p[i, j] = floatToString(p[i, j]).replace(" ", "") + "  "
                if(len(p[i, j])>mx):
                    mx = len(p[i, j])

        for i in range(p.shape[0]):
            for j in range(p.shape[1]):
                print(p[i, j] + " "*(mx-len(p[i, j])), end = "")
            print("")

        # print(p)
        properties = self.getProperties(self.last_rho, bounds)
        for prop in properties:
            if(len(prop) >3):
                print(prop[0] + " : " + floatToString(prop[1]) + " +/- " + floatToString(prop[2]))
            else:
                print(prop[0] + " : " + floatToString(prop[1]))

    # todo: this comment block
    def exportToEval(self, filePath="pythonEval.txt"):
        TORREPLACE = ""
        # Conf settings
        for k in self.conf.keys():
            if (k != "IntensityMap"):
                if (isinstance(self.conf[k], np.ndarray)):
                    A = self.conf[k]
                    TORREPLACE += "conf['" + str(k) + "'] = ["
                    for i in range(A.shape[0]):
                        TORREPLACE += str(A[i]).replace(" ", ",") + ","

                    TORREPLACE = TORREPLACE[:-1] + "]\n"
                else:
                    TORREPLACE += "conf['" + str(k) + "'] = " + str(self.conf[k]) + "\n"

        # Tomoinput
        A = self.last_input
        TORREPLACE += "tomo_input = np.array([\n"
        for i in range(A.shape[0]):
            TORREPLACE += "["
            for j in range(A.shape[1]):
                TORREPLACE += str(A[i, j]) + ","
            TORREPLACE = TORREPLACE[:-1] + "],\n"

        TORREPLACE = TORREPLACE[:-2] + "])\n"

        # intensity
        A = self.conf["IntensityMap"]
        TORREPLACE += "intensity = np.array(["
        for i in range(A.shape[0]):
            TORREPLACE += str(A[i]) + ","
        TORREPLACE = TORREPLACE[:-1] + "])"

        # print contents to file
        with open(filePath, 'w') as f:
            f.write(TORREPLACE)

    # todo: this comment block
    def exportToConf(self, filePath="pythonConf.txt"):
        TORREPLACE = ""
        # Conf settings
        for k in self.conf.keys():
            if (k != "IntensityMap"):
                if (isinstance(self.conf[k], np.ndarray)):
                    A = self.conf[k]
                    TORREPLACE += "conf['" + str(k) + "'] = ["
                    for i in range(A.shape[0]):
                        TORREPLACE += str(A[i]).replace(" ", ",") + ","

                    TORREPLACE = TORREPLACE[:-1] + "]\n"
                else:
                    TORREPLACE += "conf['" + str(k) + "'] = " + str(self.conf[k]) + "\n"

            # print contents to file
            with open(filePath, 'w') as f:
                f.write(TORREPLACE)

    # todo: this comment block
    def exportToData(self, filePath="pythonData.txt"):
        TORREPLACE = ""

        # Tomoinput
        A = self.last_input
        TORREPLACE += "tomo_input = np.array([\n"
        for i in range(A.shape[0]):
            TORREPLACE += "["
            for j in range(A.shape[1]):
                TORREPLACE += str(A[i, j]) + ","
            TORREPLACE = TORREPLACE[:-1] + "],\n"

        TORREPLACE = TORREPLACE[:-2] + "])\n"

        # intensity
        A = self.conf["IntensityMap"]
        TORREPLACE += "intensity = np.array(["
        for i in range(A.shape[0]):
            TORREPLACE += str(A[i]) + ","
        TORREPLACE = TORREPLACE[:-1] + "])"

        # print contents to file
        with open(filePath, 'w') as f:
            f.write(TORREPLACE)

    # todo: this comment block
    def exportToConf_web(self, filePath="Config_web.txt"):
        TORREPLACE = ""
        # Conf settings
        for k in self.conf.keys():
            if (k != "IntensityMap"):
                if isinstance(self.conf[k], np.ndarray):
                    A = self.conf[k]
                    TORREPLACE += "conf." + str(k) + "=["
                    for i in range(A.shape[0]):
                        TORREPLACE += str(A[i]).replace(" ", ",") + ","
                    TORREPLACE = TORREPLACE[:-1] + "];\n"
                elif isinstance(self.conf[k], list):
                    A = self.conf[k]
                    TORREPLACE += "conf." + str(k) + "=["
                    for i in range(len(A)):
                        TORREPLACE += str(A[i]).replace(" ", "") + ","
                    TORREPLACE = TORREPLACE[:-1] + "];\n"
                else:
                    TORREPLACE += "conf." + str(k) + "=" + str(self.conf[k]) + ";\n"

            # print contents to file
            with open(filePath, 'w') as f:
                f.write(TORREPLACE)

    # todo: this comment block
    def exportToData_web(self, filePath="pythonData.txt"):
        TORREPLACE = ""

        # Tomoinput
        A = self.last_input.astype("O")
        n_qubit = self.conf['NQubits']
        if self.conf['NDetectors'] == 2:
            # tomo_input[ :, 1 : 2 * n_qubit + 1 ]: singles
            A[:, 1: 2 * n_qubit + 1] = self.last_input[:, 1: 2 * n_qubit + 1].real.astype(int)
            # tomo_input[ :, 2 * n_qubit + 1 : 2 ** n_qubit + 2 * n_qubit + 1 ]: coincidences
            A[:, 2 * n_qubit + 1: 2 ** n_qubit + 2 * n_qubit + 1] = self.last_input[:,
                                                                    2 * n_qubit + 1: 2 ** n_qubit + 2 * n_qubit + 1].real.astype(
                int)
        else:
            # tomo_input[:, 1: n_qubit + 1]: singles
            A[:, 1: n_qubit + 1] = self.last_input[:, 1: n_qubit + 1].real.astype(int)
            # tomo_input[:, n_qubit + 1]: coincidences
            A[:, n_qubit + 1] = self.last_input[:, n_qubit + 1].real.astype(int)

        TORREPLACE += "tomo_input=["
        for i in range(A.shape[0]):
            TORREPLACE += "["
            for j in range(A.shape[1]):
                TORREPLACE += floatToString(A[i, j]) + ","
            TORREPLACE = TORREPLACE[:-1] + "],"

        TORREPLACE = TORREPLACE[:-1] + "];\n"

        # intensity
        A = self.conf["IntensityMap"]
        TORREPLACE += "intensity=["
        for i in range(A.shape[0]):
            TORREPLACE += str(A[i]) + ","
        TORREPLACE = TORREPLACE[:-1] + "];"

        # print contents to file
        with open(filePath, 'w') as f:
            f.write(TORREPLACE)
