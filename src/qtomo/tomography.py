import numpy as np
from . import state_utilities
import scipy

class Tomography:
    def __init__(self,
                n_qubits=2,
                n_detectors=1):
        
        self.n_qubits = n_qubits
        self.n_detectors = n_detectors

    """
    check_settings()
    desc : Checks the current settings and throws errors if any are invalid.
    """
    def check_settings(self):
        # Check Accidentals and window
        if self.n_qubits == 1:
            if self.accidental_correction:
                raise ValueError('Invalid Conf settings. Accidental Correction can not be done for single qubit tomography')
        if self.accidental_correction:
            win = np.array(self.window,dtype=float)
            if len(win.shape) != 1 or len(win) != self.n_detectors*self.n_qubits:
                raise ValueError('Invalid Conf settings. Window should have length ' +str(self.n_detectors*self.n_qubits) + " with the given settings.")
            self.window = win
        # Efficicieny and Ndetectors
        if self.n_detectors not in [1,2]:
            raise ValueError('Invalid Conf settings. NDetectors can be either 1 or 2, corresponding to the number of detectors per qubit.')
        elif self.n_detectors == 2:
            eff = np.array(self.efficiency,dtype=float)
            if not self.efficiency:
                eff = np.ones(self.n_qubits*self.n_detectors)
            if eff.ndims != 1 and eff.shape[0] != self.n_qubits*self.n_detectors:
                raise ValueError('Invalid Conf settings. Efficiency should have length ' +str(self.n_detectors*self.n_qubits) + " with the given settings.")
        elif self.n_detectors == 1:
            eff = np.ones(1)
        self.efficiency = eff
        # Crosstalk
        correctSize = int(np.floor(2**self.n_qubits+.01))
        c = self.crosstalk
        if c:
            c = np.eye(correctSize)
        else:
            c = np.array(c)
        # make sure it has the right shape
        if len(c.shape) != 2 or c.shape[0] != c.shape[1]:
            raise ValueError('Invalid Conf settings. Crosstalk should be an array with shape (' +
                                str(correctSize)+","+str(correctSize)+ ") with the given settings.")
        # make sure it has the right size
        if c.shape[0] != correctSize:
            if np.any(c-np.eye(c.shape[0])>10**-6):
                # Throw error if not the right size
                raise ValueError('Invalid Conf settings. Crosstalk should be an array with shape (' +
                                    str(correctSize) + "," + str(correctSize) + ") with the given settings.")
            else:
                # If the given cross talk is an identity matrix just fix the size.
                c = np.eye(correctSize)
        self.crosstalk = c
        

    def check_tomo_input(self, measurements, counts, intensities, rho_start, drift_correction, accidental_correction, crosstalk, efficiency, times, singles, window):

        # Check input measurements
        if measurements.ndim != 3:
            raise ValueError("Invalid measurements matrix, need an array of density matrices.")
        if measurements.shape[1] != self.n_qubits**2:
            raise ValueError("Invalid measurements matrix, the number of qubits in measurements does not match the number of qubits in the system (defined by Tomography.n_qubits).")
        if measurements.shape[0] != counts.shape[0]:
            raise ValueError("Number of counts does not match the number of measurements.")
        
        self.measurements = measurements

        # Check counts
        if (counts.shape[0] != measurements.shape[0]):
            raise ValueError("Number of Counts does not match the number of measurements")
        else:
            self.counts = counts
        # Determine if 2det/qubit from counts matrix
        if(counts.shape[1] == 1):
            self.n_detectors = 1
        elif(counts.shape[1] == self.conf['NQubits']*2):
            self.n_detectors = 2
        else:
            raise ValueError("The second axis of counts does not have the right dimension. Should be 1 or 2*NQubits for 2det")

        # Create intensities matrix
        if intensities is None:
            self.intensities = np.ones(measurements.shape[0])
        elif intensities.shape[0] == measurements.shape[0]:
            self.intensities = intensities
        else:
            raise ValueError("Invalid intensities array")
        
        # Check if drift correction is necessary
        if drift_correction or any(self.intensities-np.ones_like(self.intensities))>10**-6:
            self.drift_correction = True

        if (times or singles or window):
            self.accidental_correction = True
        else:
            self.time = np.ones(measurements.shape[0])
            self.singles = np.zeros((measurements.shape[0],self.n_detectors*self.n_qubits))
            self.window = np.zeros(self.n_detectors*self.n_qubits)


        # Check if time has right length
        if not times:
            self.time = np.ones(measurements.shape[0])
        elif not times.shape[0] == measurements.shape[0]:
            ValueError("Invalid time array")
        else:
            self.times = times

        # Check if singles has right dimensions
        if not singles:
            if self.n_detectors == 1:
                singles = np.zeros((measurements.shape[0], self.n_qubits))
            elif self.n_detectors == 2:
                singles = np.zeros((measurements.shape[0], 2 * self.n_qubits))
        else:
            if self.n_detectors == 1 and self.n_qubits == 1:
                if singles.shape != (measurements.shape[0],) and singles.shape != (measurements.shape[0], 1):
                    raise ValueError("Invalid singles matrix")
                # if the singles vector has the form (x,), this changes it to (x,1)
                if np.ndim(singles) == 1:
                    singles = np.atleast_2d(singles).transpose()
            elif self.n_detectors == 1 and self.n_qubits != 1:
                if singles.shape != (measurements.shape[0], self.n_qubits):
                    raise ValueError("Invalid singles matrix")
            elif self.n_detectors == 2:
                if singles.shape != (measurements.shape[0], 2 * self.n_qubits):
                    raise ValueError("Invalid Singles matrix")
            else:
                self.singles = singles
        # Check if window has right length
        if not window:
            self.window = np.zeros(self.n_detectors*self.n_qubits)
        elif window.n_dims == 1 and window.shape[0] != self.n_detectors*self.n_qubits:
            raise ValueError("Invalid window array")
        else:
            self.window = window

        if not crosstalk:
            self.crosstalk = np.identity(2 ** self.n_qubits)
        elif crosstalk.ndims != 2 or crosstalk.shape[0] != crosstalk.shape[1] or crosstalk.shape[0] != 2 ** self.n_qubits:
            raise ValueError("Invalid crosstalk matrix")
        else:
            self.crosstalk =  crosstalk

        if(efficiency and efficiency.shape[0] != self.n_qubits*2):
            raise ValueError("Invalid efficiency array. Length should be n_qubits*2")
        else:
            self.efficiency = efficiency

        # Initialize 
        self.rho_start = rho_start

        self.accidentals = np.zeros_like(counts)

        if intensities and efficiency:
            self.overall_norms = np.kron(intensities, efficiency)
        else:
            self.overall_norms = None

        if accidental_correction:
            scalarIndex = np.concatenate((np.ones(self.counts - 2), [2, 2]))
            additiveIndex = np.array([0, 1])
            for j in range(2, self.n_qubits):
                additiveIndex = np.concatenate(([2*j], additiveIndex))
            if (len(self.counts.shape) == 1):
                self.accidentals = self.accidentals[:, np.newaxis]
            for j in range(self.counts.shape[0]):
                index = bin(j).split("b")[1]
                index = "0" * (self.n_qubits_ - len(index)) + index
                index = [int(char) for char in index]
                index = index*scalarIndex + additiveIndex
                index = np.array(index, dtype = int)
                self.accidentals[:, j] = np.prod(np.real(self.singles[:, tuple(index)]), axis = 1) * (window[j] * 1e-9 / np.real(self.times)) ** (self.n_qubits - 1)
            if (self.accidentals.shape != self.counts.shape):
                self.accidentals = self.accidentals[:, 0]
    

    def tomography(self, measurements, counts, rho_start = None, crosstalk=None, accidental_correction=False, drift_correction = False, efficiency=None, times=None, singles=None, window=None, error=0,
                        intensities=None, method="MLE"):
        
        self.check_tomo_input(measurements, counts, intensities, rho_start, drift_correction, accidental_correction, crosstalk, efficiency, times, singles, window)
        self.check_settings()
        # get the starting state from tomography_LINEAR if not defined
        if self.rho_start is None:
            try:
                [self.rho_start,inten_linear] = self.tomography_LINEAR(self.counts, self.measurements, self.overall_norms)
                # Currently linear tomography gets the phase wrong. So a temporary fix is to just transpose it.
                self.rho_start = self.rho_start.transpose()
                self.rho_start = state_utilities.make_positive(self.rho_start)
                self.rho_start = self.rho_start / np.trace(self.rho_start)
            except Exception:
                raise RuntimeError('Failed to run linear Tomography', str(Exception))

        # Run tomography and find an estimate for the state
        if method == "linear":
            [self.rhog, self.intensity, self.fvalp] = [self.rho_start,inten_linear,0]
        elif method == "mle":
            [self.rhog, self.intensity, self.fvalp] = self.tomography_MLE(self.rho_start, self.counts, self.measurements, self.accidentals, self.overall_norms)
        elif method == "hmle":
            [self.rhog, self.intensity, self.fvalp] = self.tomography_HMLE(self.rho_start, self.counts, self.measurements, self.accidentals, self.overall_norms)
        # elif method.upper() == "BME":
        #     [rhog, intensity, fvalp] = self.tomography_BME(starting_matrix, coincidences, measurements_densities,accidentals,overall_norms)
        else:
            raise ValueError("Invalid Method name: " + str(method))

        self.mont_carlo_states = list([[self.rhog, self.intensity, self.fvalp]])

        return [self.rhog, self.intensity, self.fvalp]
        
    """
    tomography_LINEAR(coincidences, measurements,overall_norms)
    Desc: Uses linear techniques to find a starting state for maximum likelihood estimation.

    Parameters
    ----------
    coincidences : ndarray shape = (Number of measurements, NDetectors**NQubits)
        The counts of the tomography.
    measurements : ndarray with shape = (NDetectors*number of measurements,2^numQubits) 
        The measurements of the tomography in pure state form.
    overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The relative weights of each measurment. Used for drift correction.

    Returns
    -------
    rhog : ndarray with shape = (2^numQubits, 2^numQubits)
        The starting predicted state.
    intensity : float
        The predicted overall intensity used to normalize the state.
    """

    def tomography_LINEAR(self, counts, measurements, overall_norms=None):
        # If overall_norms not given then assume uniform
        if overall_norms is None:
            overall_norms = np.ones(counts.shape[0])
        elif not overall_norms.shape[0] == counts.shape[0]:
            raise ValueError("Invalid intensities array")

        counts = counts.flatten()

        pauli_basis = state_utilities.generalized_pauli_basis(self.getNumQubits())
        stokes_measurements = np.array([state_utilities.get_stokes_parameters(m,pauli_basis) for m in measurements]) / 2**self.getNumQubits()
        freq_array = counts / overall_norms

        B_inv = np.linalg.inv(np.matmul(stokes_measurements.T,stokes_measurements))
        stokes_params = np.matmul(stokes_measurements.T,freq_array)
        stokes_params = np.matmul(B_inv,stokes_params)
        linear_rhog = np.multiply(pauli_basis,stokes_params[:, np.newaxis,np.newaxis])
        linear_rhog = np.sum(linear_rhog,axis=0)

        intensity = np.trace(linear_rhog)
        rhog = linear_rhog / intensity

        return [rhog, intensity]
    
    """
    tomography_MLE(starting_matrix, coincidences, measurements, accidentals,overall_norms)
    Desc: Runs tomography using maximum likelyhood estimation.

    Parameters
    ----------
    starting_matrix : ndarray with shape = (2^numQubits, 2^numQubits)
        The starting predicted state.
    coincidences : ndarray shape = (Number of measurements, NDetectors**NQubits)
        The counts of the tomography.
    measurements_densities : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
        The measurements of the tomography in density matrix form.
    accidentals : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The accidental values of the tomography. Used for accidental correction.
    overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The relative weights of each measurment. Used for drift correction.

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
    def tomography_MLE(self, starting_matrix, coincidences, measurements, accidentals,overall_norms=-1):
        # If overall_norms not given then assume uniform
        if not isinstance(overall_norms,np.ndarray):
            overall_norms = np.ones(coincidences.shape[0])
        elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
            raise ValueError("Invalid intensities array")

        init_intensity = np.mean(np.multiply(coincidences, 1/overall_norms)) * (starting_matrix.shape[0])
        starting_tvals = state_utilities.density2t(starting_matrix)
        starting_tvals = starting_tvals + 0.0001
        starting_tvals = starting_tvals * np.sqrt(init_intensity)

        coincidences = np.real(coincidences)
        coincidences = coincidences.flatten()

        final_tvals = scipy.optimize.leastsq(maxlike_fitness, np.real(starting_tvals),
                              args = (coincidences, accidentals, measurements, overall_norms),
                              ftol=self.conf["ftol"],
                              xtol=self.conf["xtol"],
                              gtol=self.conf["gtol"],
                              maxfev=self.conf["maxfev"])[0]
        fvalp = np.sum(maxlike_fitness(final_tvals, coincidences, accidentals, measurements, overall_norms) ** 2)

        final_matrix = state_utilities.t_to_density(final_tvals, normalize=False)
        intensity = np.trace(final_matrix)
        final_matrix = final_matrix / np.trace(final_matrix)

        intensity = np.float64(np.real(intensity))

        return [final_matrix, intensity, fvalp]
    
    """
    tomography_HMLE(starting_matrix, coincidences, measurements, accidentals,overall_norms)
    Desc: Runs tomography using hedged maximum likelyhood estimation.

    Parameters
    ----------
    starting_matrix : ndarray with shape = (2^numQubits, 2^numQubits)
        The starting predicted state.
    coincidences : ndarray shape = (Number of measurements, NDetectors**NQubits)
        The counts of the tomography.
    measurements_densities : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
        The measurements of the tomography in density matrix form.
    accidentals : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The accidental values of the tomography. Used for accidental correction.
    overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The relative weights of each measurment. Used for drift correction.

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
    def tomography_HMLE(self, starting_matrix, coincidences, measurements, accidentals,overall_norms=-1):
        # If overall_norms not given then assume uniform
        if not isinstance(overall_norms, np.ndarray):
            overall_norms = np.ones(coincidences.shape[0])
        elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
            raise ValueError("Invalid intensities array")

        init_intensity = np.mean(np.multiply(coincidences, 1 / overall_norms)) * (starting_matrix.shape[0])
        starting_tvals = state_utilities.density2t(starting_matrix)
        starting_tvals = starting_tvals + 0.0001
        starting_tvals = starting_tvals * np.sqrt(init_intensity)

        coincidences = np.real(coincidences)
        coincidences = coincidences.flatten()

        bet = self.conf['Beta']
        if bet > 0:
            final_tvals = scipy.optimize.leastsq(maxlike_fitness_hedged, np.real(starting_tvals),
                                  args=(coincidences, accidentals, measurements, bet, overall_norms),
                                  ftol=self.conf["ftol"],
                                  xtol=self.conf["xtol"],
                                  gtol=self.conf["gtol"],
                                  maxfev=self.conf["maxfev"])[0]
            fvalp = np.sum(maxlike_fitness_hedged(final_tvals, coincidences, accidentals, measurements, bet,
                                                      overall_norms) ** 2)
        else:
            raise ValueError("To use Hedged Maximum Likelihood, Beta must be a positive number.")

        final_matrix = state_utilities.t_to_density(final_tvals,normalize=False)
        intensity = np.trace(final_matrix)
        final_matrix = final_matrix / np.trace(final_matrix)

        intensity = np.float64(np.real(intensity))

        return [final_matrix, intensity, fvalp]
    
    '''
    tomography_BME(starting_matrix, coincidences, measurements, accidentals,overall_norms)
    Desc: Runs tomography using bayesian mean estimation. Currently in progress.

    Parameters
    ----------
    starting_matrix : ndarray with shape = (2^numQubits, 2^numQubits)
        The starting predicted state.
    coincidences : ndarray shape = (Number of measurements, NDetectors**NQubits)
        The counts of the tomography.
    measurements : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
        The measurements of the tomography in density matrix form.
    accidentals : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The accidental values of the tomography. Used for accidental correction.
    overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The relative weights of each measurment. Used for drift correction.

    Returns
    -------
    rhog : ndarray with shape = (2^numQubits, 2^numQubits)
        The predicted density matrix.
    intensity : float
        The predicted overall intensity used to normalize the state.
    fvalp : float
        Final value of the internal optimization function. Values greater than the number
        of measurements indicate poor agreement with a quantum state.
    '''
    # def tomography_BME(self,starting_matrix, coincidences, measurements, accidentals,overall_norms):
    #
    #     # algorithm parameters
    #     preallocationSize = 50000  # preallocate memory for speed
    #     changeThreshold = 10 ** -19 # stopping condition threshold
    #     numberOfPriorSamples = 10000  # min number of samples to sample from the inital prior
    #     updateFrequency = 300  # controls how often the estimate of the posterior is updated, and how often stability is checked
    #     max_iterations = 98  # The maximum number of iterations(updates)
    #
    #     # initialize
    #     d = starting_matrix.shape[0]
    #     i = 0
    #     iterationCounter = 0
    #     eps = 2.2204e-16
    #     prevMean_rho = np.zeros((d,d))
    #     prev_minIndex = -1
    #     mean_tvals = -1
    #     cov_tvals = -1
    #
    #     # Todo TEMP
    #     tempLogLikeMean = np.zeros(101)
    #     tempLogLikeMax = np.zeros(101)
    #     tempLogLikeSD = np.zeros(101)
    #
    #
    #     # Lists that are evaluated periodically
    #     stabilityHistory = np.zeros(max_iterations*3)
    #
    #     # List of Sample States
    #     sampleStates_rho = np.zeros((d,d,preallocationSize), dtype=complex)
    #     sampleStates_loglike = np.zeros(preallocationSize)
    #     sampleStates_tvals = np.zeros((preallocationSize, d ** 2))
    #
    #     # estimate the counts intensity
    #     norm = sum(coincidences / (len(coincidences) * d))
    #     minimization_output = minimize(log_likelyhood, norm,
    #                           args=(density2t(starting_matrix), coincidences, measurements, accidentals,overall_norms))
    #     norm = minimization_output.x[0]
    #     base_loglike = minimization_output.fun
    #
    #     # SAMPLE FROM PRIOR
    #     # Samples are drawn from the ginibre distribution
    #     stillUsingPrior = 1
    #     while stillUsingPrior:
    #         # expand memory in chunks if preallocation size is exceeded
    #         if (i >= len(sampleStates_loglike)):
    #             sampleStates_rho = np.concatenate((sampleStates_rho, np.zeros((d,d,i))), axis=2)
    #             sampleStates_loglike = np.concatenate((sampleStates_loglike, np.zeros(i)))
    #             sampleStates_tvals = np.concatenate((sampleStates_tvals, np.zeros((i, d ** 2))), axis=0)
    #
    #         # Sample from the ginibre distribution
    #         random_rho = random_density_state(self.conf['NQubits'])
    #         random_tvals = density2t(random_rho)
    #         random_loglike = log_likelyhood(norm, random_rho, coincidences, measurements, accidentals,overall_norms)
    #
    #         # estimate the counts intensity
    #         if random_loglike < base_loglike:
    #             minimization_output = minimize(log_likelyhood, norm,
    #                                            args=(random_rho, coincidences, measurements,
    #                                                  accidentals,
    #                                                  overall_norms))
    #             norm = minimization_output.x[0]
    #             random_loglike = minimization_output.fun
    #             base_loglike = random_loglike
    #
    #         # Store sample state, tval, and its loglikelyhood
    #         sampleStates_rho[:,:,i] = random_rho
    #         sampleStates_tvals[i] = random_tvals
    #         sampleStates_loglike[i] = random_loglike
    #
    #
    #         # Periodically perform the following tasks
    #         if i % updateFrequency == 0 and i !=0:
    #             # Normalize the likelihoods
    #             [sampleStates_normlike, minIndex,sampleStates_scaledloglike] = normalizeLikelihoods(sampleStates_loglike[:i + 1])
    #             unNormedMean_rho = np.dot(sampleStates_rho[:, :, :i + 1], np.exp(-1*sampleStates_scaledloglike))
    #
    #             # Update the mean
    #             tr = np.trace(unNormedMean_rho)
    #             if tr == 0:
    #                 mean_rho = np.zeros_like(prevMean_rho)
    #                 stabilityHistory[iterationCounter] = float('inf')
    #             else:
    #                 mean_rho = unNormedMean_rho / tr
    #                 if iterationCounter != 0:
    #                     # Calculate the stability
    #                     sectionLikelihood = sum(np.exp(-1 * sampleStates_scaledloglike[(i + 1 - updateFrequency):(i + 1)]))
    #                     totalLikelihood = sum(np.exp(-1 * sampleStates_scaledloglike[1:(i + 1)]))
    #                     normFactor = totalLikelihood/((iterationCounter+1)*sectionLikelihood)
    #                     infidelity = (1 - np.real(fidelity(mean_rho, prevMean_rho)))
    #                     # ensure rounding never leads the infidelity to be negative
    #                     if infidelity < 0:
    #                         infidelity = eps
    #                     stabilityHistory[iterationCounter] = infidelity * normFactor
    #                 else:
    #                     stabilityHistory[iterationCounter] = 0
    #             prevMean_rho = mean_rho
    #
    #             # If we have met the min required amount of prior samples try to switching to the posterior estimate
    #             if i >= numberOfPriorSamples:
    #                 contributingParams = sampleStates_normlike > 0
    #                 firstMax_like = max(sampleStates_normlike)
    #                 secondMax_like = max(sampleStates_normlike[sampleStates_normlike < firstMax_like])
    #
    #                 # IMPOSE STOPPING CONDITION IF TOO MANY ITERATIONS REACHED
    #                 if (iterationCounter > max_iterations):
    #                     warnings.warn(
    #                         "MAX ITER REACHED - BME (Prior). Stability: " + str(
    #                             stabilityHistory[iterationCounter]))
    #                     try:
    #                         [mean_tvals, cov_tvals] = weightedcov(sampleStates_tvals[:i + 1][contributingParams, :],
    #                                                               sampleStates_normlike[contributingParams])
    #                     except:
    #                         pass
    #                     stillUsingPrior = 0
    #                     i += 1
    #                     iterationCounter = iterationCounter + 1
    #                     break
    #
    #                 # Switch to the posterior under the right conditions
    #                 # 1.) We must ensure one of our samples does not make up majority of the contribution in estimating
    #                 #     the posterior parameters
    #                 # 2.) Cov must be positive semi - definite
    #                 if (secondMax_like / firstMax_like) >= .01:
    #                     [mean_tvals, cov_tvals] = weightedcov(sampleStates_tvals[:i + 1][contributingParams, :],
    #                                                           sampleStates_normlike[contributingParams])
    #                     # Try to find eigenvals. If it can't then keep sampling.
    #                     try:
    #                         if all(np.linalg.eigvals(cov_tvals) >= 0):
    #                             # Stop using the prior
    #                             stillUsingPrior = 0
    #                             i += 1
    #                             iterationCounter = iterationCounter + 1
    #                             break
    #                     except:
    #                         pass
    #             # Increase iteration number
    #             iterationCounter = iterationCounter + 1
    #         i = i + 1
    #
    #
    #     iterPrior = iterationCounter
    #     self.priorsToConverge = i-1
    #
    #
    #     # todo TEMP
    #     tempLogLikeMean[0] = np.average(normalizeLikelihoods(sampleStates_loglike[:i])[0])
    #     tempLogLikeMax[0] = np.max(normalizeLikelihoods(sampleStates_loglike[:i])[0])
    #     tempLogLikeSD[0] = np.std(normalizeLikelihoods(sampleStates_loglike[:i])[0])
    #
    #
    #     # SAMPLE FROM POSTERIOR
    #     # Samples by drawing tVals from a multivariate normal distro
    #     stillUsingPosterior = 1
    #     while stillUsingPosterior:
    #         # expand memory in chunks if preallocation size is exceeded
    #         if (i >= len(sampleStates_loglike)):
    #             sampleStates_rho = np.concatenate((sampleStates_rho, np.zeros((d,d,i))), axis=2)
    #             sampleStates_loglike = np.concatenate((sampleStates_loglike, np.zeros(i)))
    #             sampleStates_tvals = np.concatenate((sampleStates_tvals, np.zeros((i, d ** 2))), axis=0)
    #
    #         # Draws by sampling the tVals from a multivariate normal distro
    #         random_tvals = rand.multivariate_normal(mean_tvals, cov_tvals)
    #         random_rho = t_to_density(random_tvals)
    #         random_loglike = log_likelyhood(norm, random_rho, coincidences, measurements, accidentals,overall_norms)
    #
    #         # estimate the counts intensity
    #         if random_loglike < base_loglike:
    #             minimization_output = minimize(log_likelyhood, norm,
    #                                            args=(random_rho, coincidences, measurements,
    #                                                  accidentals,
    #                                                  overall_norms))
    #             norm = minimization_output.x[0]
    #             random_loglike = minimization_output.fun
    #             base_loglike = random_loglike
    #
    #         # Store sample state, tval, and its loglikelyhood
    #         sampleStates_rho[:, :, i] = random_rho
    #         sampleStates_tvals[i] = random_tvals
    #         sampleStates_loglike[i] = random_loglike
    #
    #         # Periodically perform the following tasks
    #         if i % updateFrequency == 0:
    #             # Normalize the likelihoods
    #             [sampleStates_normlike, minIndex, sampleStates_scaledloglike] = normalizeLikelihoods(
    #                 sampleStates_loglike[:i + 1])
    #             unNormedMean_rho = np.dot(sampleStates_rho[:, :, :i + 1], np.exp(-1*sampleStates_scaledloglike))
    #
    #
    #             # todo TEMP
    #             tempLogLikeMean[iterationCounter-iterPrior+1] = np.average(sampleStates_normlike)
    #             tempLogLikeMax[iterationCounter-iterPrior+1] = max(sampleStates_normlike)
    #             tempLogLikeSD[iterationCounter-iterPrior+1] = np.std(sampleStates_normlike)
    #
    #
    #
    #             # Update the mean
    #             tr = np.trace(unNormedMean_rho)
    #             if tr == 0:
    #                 mean_rho = np.zeros_like(prevMean_rho)
    #                 stabilityHistory[iterationCounter] = float('inf')
    #             else:
    #                 mean_rho = unNormedMean_rho / tr
    #
    #                 # Calculate the stability
    #                 sectionLikelihood = sum(np.exp(-1*sampleStates_scaledloglike[(i + 1 - updateFrequency):(i + 1)]))
    #                 totalLikelihood = sum(np.exp(-1*sampleStates_scaledloglike[1:(i + 1)]))
    #                 normFactor = totalLikelihood / ((iterationCounter + 1) * sectionLikelihood)
    #                 infidelity = (1 - np.real(fidelity(mean_rho, prevMean_rho)))
    #                 # ensure rounding never leads the infidelity to be negative
    #                 if infidelity < 0:
    #                     infidelity = eps
    #                 stabilityHistory[iterationCounter] = infidelity * normFactor
    #
    #                 # Check if the stopping condition is met
    #                 if (iterationCounter > iterPrior) and (stabilityHistory[iterationCounter] < changeThreshold):
    #                     stillUsingPosterior = 0
    #
    #             prevMean_rho = mean_rho
    #
    #             # IMPOSE STOPPING CONDITION IF TOO MANY ITERATIONS REACHED
    #             if (iterationCounter > (max_iterations+iterPrior)):
    #                 warnings.warn(
    #                     "MAX ITER REACHED - BME (Post). Stability: " + str(stabilityHistory[iterationCounter]))
    #                 stillUsingPosterior = 0
    #
    #             contributingParams = sampleStates_normlike > 0
    #             firstMax_like = max(sampleStates_normlike)
    #             secondMax_like = max(sampleStates_normlike[sampleStates_normlike < firstMax_like])
    #
    #             # Update the posterior estimate only if:
    #             # 1.) We must ensure one of our samples make up majority of the contribution in estimating
    #             #     the posterior parameters
    #             # 2.) Cov must be positive semi - definite
    #             if (secondMax_like / firstMax_like) > .001:
    #                 [newMean_tvals, newCov_tvals] = weightedcov(sampleStates_tvals[:i + 1][contributingParams, :],
    #                                                             sampleStates_normlike[contributingParams])
    #                 try:
    #                     # Try to find eigenvals. If it can't then keep sampling.
    #                     if all(np.linalg.eigvals(newCov_tvals) >= 0):
    #                         cov_tvals = newCov_tvals
    #                         mean_tvals = newMean_tvals
    #                 except:
    #                     pass
    #             # Increase iteration number
    #             iterationCounter = iterationCounter + 1
    #         # Increase sample number
    #         i = i + 1
    #     self.stabilityHistory = stabilityHistory
    #     print("done.",end="")
    #     samplesToSC = i
    #
    #
    #     self.tempData = [tempLogLikeMean,tempLogLikeSD,tempLogLikeMax]
    #     return [mean_rho, norm, samplesToSC]

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
    intenp : 1darray with length = n
        The intensity of each approximate density matrices.
    fvalp : 1darray with length = n
        The fval of the regression associated with each approximate density matrices.
    """
    def tomography_states_generator(self, n, method="mle"):
        # Save the last data so we can restore it later
        last_input = self.last_input.copy()

        acc = self.accidental_correction
        time = self.times
        meas = self.measurements
        singles = self.singles
        counts = self.counts

        # Re-sample counts and singles
        test_counts = np.random.poisson(counts)
        test_singles = np.random.poisson(singles)
        for j in range(n):
            if len(test_counts.shape) == 1:
                test_counts = np.array([test_counts]).T
            [rhop, intenp, fvalp] = self.tomography(meas,test_counts,times=time, accidental_correction = acc, singles=test_singles,intensities=self.intensities,method=method)
            self.mont_carlo_states.append([rhop, intenp, fvalp])
        # Restore the last tomo_input matrix to the original one
        self.last_input = last_input
        return self.mont_carlo_states
    
    """
    getBellSettings(bounds)
    Desc: Returns the optimal measurment settings for the CHSH bell inequality. In Progress, have not checked.
    """
    def getBellSettings(self, partsize_init = 9, partsize = 5, t = 3, bounds = -1):
        if self.conf['nqubits'] != 2:
            raise ValueError('Invalid Conf settings. Bell state can only be found for 2 qubit tomography')
        rho = self.last_rho
        # if bounds not set use the conf settings
        if bounds == -1:
            bounds = self.conf['DoErrorEstimation']

        # generate states if needed
        if bounds > len(self.mont_carlo_states) - 1:
            self.tomography_states_generator(bounds - len(self.mont_carlo_states) + 1)
        states = np.array(self.mont_carlo_states,dtype="O")[:bounds + 1, :]

        if (bounds > 0):
            return getBellSettings_helper_bounds(states[:,0], rho, partsize_init, partsize, t, bounds)
        else:
            return getBellSettings_helper(rho, partsize_init, partsize, t)

# Helper function for calculating the bell settings
def coinmat(a, b):
    k = np.array(
        [
            np.cos(a) * np.cos(b),
            np.cos(a) * np.sin(b),
            np.sin(a) * np.cos(b),
            np.sin(a) * np.sin(b),
        ]
    )
    cmat = np.outer(k, k)

    return cmat

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
overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
    The relative weights of each measurment. Used for drift correction.
    
Returns
-------
val : float
    value of the optimization function.
"""


def maxlike_fitness(t, coincidences, accidentals, measurements, overall_norms):
    rhog = state_utilities.t_to_density(t, normalize=False)
    prediction = np.zeros_like(coincidences)
    for j in range(len(prediction)):
        prediction[j] = (
            overall_norms[j] * np.real(np.trace(np.dot(measurements[j, :, :], rhog)))
            + accidentals[j]
        )
        prediction[j] = np.max([prediction[j], 0.01])
    log_like = (prediction - coincidences) / np.sqrt(prediction)
    return np.real(log_like)


"""
maxlike_fitness_hedged(t, coincidences, accidentals, measurements, prediction, bet)
Desc: Calculates the diffrence between the current predicted state data and the actual data using hedged maximum likelihood.

Parameters
----------
t : ndarray
    T values of the current predicted state.
coincidences : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
    The counts of the tomography.
accidentals : ndarray with length = number of measurements or shape = (number of measurements, 2^numQubits) for 2 det/qubit
    The singles values of the tomography. Used for accidental correction.
measurements : ndarray with shape = (2^numQubits, 2^numQubits, number of measurements)
    The measurements of the tomography in density matrix form .
prediction : ndarray
    Predicted counts from the predicted state.
bet : float 0 to .5
    The beta value used.
overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
    The relative weights of each measurment. Used for drift correction.

Returns
-------
val : float
    value of the optimization function.
"""


def maxlike_fitness_hedged(
    t, coincidences, accidentals, measurements, bet, overall_norms
):
    prediction = np.zeros_like(coincidences)
    rhog = state_utilities.t_to_density(t, normalize=False)
    for j in range(len(prediction)):
        prediction[j] = (
            overall_norms[j] * np.real(np.trace(np.dot(measurements[j, :, :], rhog)))
            + accidentals[j]
        )
        prediction[j] = np.max([prediction[j], 0.01])
    hedge = np.repeat(
        np.real((bet * np.log(np.linalg.det(rhog))) / len(prediction)), len(prediction)
    )
    val = np.sqrt(
        np.real((((prediction - coincidences) ** 2) / (2 * prediction)) - hedge) + 1000
    )
    return np.real(val)


# used in getBellSettings
def bellsettings_range_init(rhog, partsize):
    sval = 0
    aval = 0
    apval = 0
    bval = 0
    bpval = 0

    for a in np.linspace(0, np.pi / 2, partsize):
        for ap in np.linspace(a, np.pi / 2, partsize):
            for b in np.linspace(0, np.pi / 2, partsize):
                for bp in np.linspace(b, np.pi / 2, partsize):
                    npp = np.real(np.trace(np.dot(coinmat(a, b), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(a + np.pi / 2, b + np.pi / 2), rhog))
                    )
                    e_ab = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, b), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(ap + np.pi / 2, b + np.pi / 2), rhog))
                    )
                    e_apb = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(a, bp), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(a + np.pi / 2, bp + np.pi / 2), rhog))
                    )
                    e_abp = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, bp), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(ap + np.pi / 2, bp + np.pi / 2), rhog))
                    )
                    e_apbp = 2 * (npp + nmm) - 1

                    s = (
                        e_ab
                        + e_abp
                        + e_apb
                        + e_apbp
                        - 2 * np.min([e_ab, e_abp, e_apb, e_apbp])
                    )

                    if s > sval:
                        sval = s
                        aval = a
                        apval = ap
                        bval = b
                        bpval = bp

    arange_s = [
        np.max([aval - ((np.pi / 2) / partsize), 0]),
        np.min([aval + ((np.pi / 2) / partsize), np.pi / 2]),
    ]
    aprange_s = [
        np.max([apval - ((np.pi / 2) / partsize), 0]),
        np.min([apval + ((np.pi / 2) / partsize), np.pi / 2]),
    ]
    brange_s = [
        np.max([bval - ((np.pi / 2) / partsize), 0]),
        np.min([bval + ((np.pi / 2) / partsize), np.pi / 2]),
    ]
    bprange_s = [
        np.max([bpval - ((np.pi / 2) / partsize), 0]),
        np.min([bpval + ((np.pi / 2) / partsize), np.pi / 2]),
    ]

    return [sval, arange_s, brange_s, aprange_s, bprange_s]


# used in getBellSettings
def bellsettings_range(rhog, partsize, arange, brange, aprange, bprange):
    sval = 0
    aval = 0
    apval = 0
    bval = 0
    bpval = 0

    for a in np.linspace(arange[0], arange[1], partsize):
        for ap in np.linspace(aprange[0], aprange[1], partsize):
            for b in np.linspace(brange[0], brange[1], partsize):
                for bp in np.linspace(bprange[0], bprange[1], partsize):
                    npp = np.real(np.trace(np.dot(coinmat(a, b), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(a + np.pi / 2, b + np.pi / 2), rhog))
                    )
                    e_ab = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, b), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(ap + np.pi / 2, b + np.pi / 2), rhog))
                    )
                    e_apb = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(a, bp), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(a + np.pi / 2, bp + np.pi / 2), rhog))
                    )
                    e_abp = 2 * (npp + nmm) - 1

                    npp = np.real(np.trace(np.dot(coinmat(ap, bp), rhog)))
                    nmm = np.real(
                        np.trace(np.dot(coinmat(ap + np.pi / 2, bp + np.pi / 2), rhog))
                    )
                    e_apbp = 2 * (npp + nmm) - 1

                    s = (
                        e_ab
                        + e_abp
                        + e_apb
                        + e_apbp
                        - 2 * np.min([e_ab, e_abp, e_apb, e_apbp])
                    )

                    if s > sval:
                        sval = s
                        aval = a
                        apval = ap
                        bval = b
                        bpval = bp

    arange_s = [
        np.max([aval - ((arange[1] - arange[0]) / partsize), 0]),
        np.min([aval + ((arange[1] - arange[0]) / partsize), np.pi / 2]),
    ]
    aprange_s = [
        np.max([apval - ((aprange[1] - aprange[0]) / partsize), 0]),
        np.min([apval + ((aprange[1] - aprange[0]) / partsize), np.pi / 2]),
    ]
    brange_s = [
        np.max([bval - ((brange[1] - brange[0]) / partsize), 0]),
        np.min([bval + ((brange[1] - brange[0]) / partsize), np.pi / 2]),
    ]
    bprange_s = [
        np.max([bpval - ((bprange[1] - bprange[0]) / partsize), 0]),
        np.min([bpval + ((bprange[1] - bprange[0]) / partsize), np.pi / 2]),
    ]

    return [sval, aval, apval, bval, bpval, arange_s, brange_s, aprange_s, bprange_s]


# Helper function that returns the optimal CHSH bell measurement settings given rho.
# It is recommended to use getBellSettings in the TomoClass
def getBellSettings_helper(rhog, partsize_init, partsize, t):
    [s, arange, brange, aprange, bprange] = bellsettings_range_init(rhog, partsize_init)
    a = 0
    b = 0
    ap = 0
    bp = 0

    for j in range(t):
        [s, a, ap, b, bp, arange, brange, aprange, bprange] = bellsettings_range(
            rhog, partsize, arange, brange, aprange, bprange
        )

    return np.array([["s", s], ["a", a], ["a'", ap], ["b", b], ["b'", bp]], dtype="O")


# Helper function that returns the optimal CHSH bell measurement settings with bounds of the given rhos.
# It is recommended to use getBellSettings in the TomoClass
def getBellSettings_helper_bounds(rhop, rho, partsize_init, partsize, t, n):
    belldata = np.zeros([n + 1, 5])

    for j in range(n):
        belldata[j] = getBellSettings_helper(rhop[j], partsize_init, partsize, t)[:, 1]
    [bellNames, belldata[-1]] = getBellSettings_helper(
        rho, partsize_init, partsize, t
    ).transpose()
    bmeans = np.zeros(5)
    berrors = np.zeros(5)

    for m in range(5):
        berrors[m] = np.std(belldata[:, m], ddof=n - 1)
        bmeans[m] = np.mean(belldata[:, m])

    return np.array([bellNames, berrors, bmeans], dtype="O").transpose()

    