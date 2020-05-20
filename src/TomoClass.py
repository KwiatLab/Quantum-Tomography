from .TomoFunctions import *
from .TomoHelpers import *
import numpy as np
from scipy.optimize import leastsq



"""
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """


"""
    class Tomography()
    Desc: This is the main tomography object that the library is built around. The goal is to only have one tomography object and edit the 
    configuration settings as you go. The conf settings are saved in the class however the input data and the output data are not stored in the object.
    """
class Tomography():

    ######################
    '''Public Variables'''
    ######################

    # conf['NQubits']: >= 1, it will take a much longer time for more qubits.
    # conf['NDetectors']: 1 or 2
    # conf['ctalk']: [[C0->0, C0->1],[C1->0, C1->1]]
    # conf['Bellstate']: 'no' or 'yes'
    # conf['DoDriftCorrection'] = 'no' or 'yes'
    # conf['DoAccidentalCorrection'] = 'no' or 'yes'
    # conf['DoErrorEstimation']: >= 0
    # conf['Window']: 0 or array like, dimension = 1
    # conf['Efficiency']: 0 or array like, dimension = 1
    # conf['Beta']: 0 to 0.5, depending on purity of state and total number of measurements.

    # tomo_input: array like, dimension = 2.
    #
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

    ##################
    '''Constructors'''
    ##################

    """
    Default Constructor
    Desc: This initializes a default tomography object
    """
    def __init__(self):
        self.conf = {'NQubits': 2,
            'NDetectors': 1,
            'Crosstalk': np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]),
            'Bellstate': 0,
            'DoDriftCorrection': 0,
            'DoAccidentalCorrection' : 1,
            'DoErrorEstimation': 0,
            'Window': 0,
            'Efficiency': 0,
            'RhoStart': [],
            'IntensityMap': [[1]],
            'Beta': 0}
        self.err_functions = ['concurrence','tangle','entanglement','entropy','linear_entropy','negativity','purity']

    """
    setConfSetting(setting,val)
    Desc: Sets a specific self.conf setting

    Parameters
    ----------
    setting : string
        The setting you want to update. 
        Possible values are ['NQubits','NDetectors','Crosstalk','Bellstate','DoDriftCorrection','DoAccidentalCorrection','DoErrorEstimation','Window','Efficiency','RhoStart','IntensityMap','Beta']
    val: array_like, int, or string
            The new value you want to the setting to be.
    """
    def setConfSetting(self,setting,val):
        try:
            valC = val.copy()
        except:
            valC = val
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
    def importConf(self,conftxt):
        conf = self.conf
        exec(compile(open(conftxt, "rb").read(), conftxt, 'exec'))
        self.standardizeConf()

    """
   importData(datatxt)
   Desc: Import a text file containing the tomography data and run tomography.

   Parameters
   ----------
   datatxt : string
       path to data file
   """
    def importData(self,datatxt):
        exec(compile(open(datatxt, "rb").read(), datatxt, 'exec'))
        return self.state_tomography(locals().get('tomo_input'), locals().get('intensity'))

    """
    importEval(evaltxt)
    Desc: Import a eval file containing the tomography data and the configuration setting, then run tomography.

    Parameters
    ----------
    evaltxt : string
        path to eval file
    """
    def importEval(self,evaltxt):
        conf = self.conf
        exec(compile(open(evaltxt, "rb").read(), evaltxt, 'exec'))
        self.standardizeConf()
        return self.state_tomography(locals().get('tomo_input'), locals().get('intensity'))
    """
    standardizeConf()
    Desc: Helper function to handle different cases of conf inputs
    """
    def standardizeConf(self):
        for k in self.conf.keys():
            if (isinstance(self.conf[k], str)):
                if (self.conf[k].lower() == "yes" or self.conf[k].lower() == "true"):
                    self.conf[k] = 1
                elif (self.conf[k].lower() == "no" or self.conf[k].lower() == "false"):
                    self.conf[k] = 0

    ##########################
    '''Tomography Functions'''
    ##########################

    """
    state_tomography(raw_counts, intensities)
    Desc: Main function that runs tomography.
    
    Parameters
    ----------
    raw_counts : ndarray
        The input data for the current tomography. This is what tomo_input will be set to.
    intensities : ndarray with length = number of measurements
        Relative pump power (arb. units) during measurement; used for drift correction.
    Returns
    -------
    rhog : ndarray with size = (2^numQubits,2^numQubits) 
        The predicted density matrix.
    intensity : float
        The predicted overall intensity used to normalize the state.
    fvalp : float
        Final value of the internal optimization function. Values greater than the number
        of measurements indicate poor agreement with a quantum state.
    """
    def state_tomography(self,raw_counts, intensities):
        rho0 = self.conf['RhoStart']
        self.tomo_input = raw_counts
        self.intensity = intensities
        [data, m1, m2, acc] = self.filter_data(raw_counts, intensities)

        if not rho0:
            rho0 = self.linear_tomography(data, m2)[0]

        [rhog, intensity, fvalp] = self.maximum_likelihood_tomography(rho0, data, m1, acc)
        return [rhog, intensity, fvalp]


    """
    maximum_likelihood_tomography(rho0, data, m, acc)
    Desc: Calculates the most likely state given the data.

    Parameters
    ----------
    rho0 : ndarray with size = (2^numQubits,2^numQubits) 
        The starting predicted state found with linear tomography.
    data : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The counts of the tomography.
    m : ndarray with size = (2^numQubits,2^numQubits,number of measurements) 
        The measurements of the tomography in density matrix form.
    acc : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.
    Returns
    -------
    rhog : ndarray with size = (2^numQubits,2^numQubits) 
        The predicted density matrix.
    intensity : float
        The predicted overall intensity used to normalize the state.
    fvalp : float
        Final value of the internal optimization function. Values greater than the number
        of measurements indicate poor agreement with a quantum state.
    """
    def maximum_likelihood_tomography(self,rho0, data, m, acc):
        rho0 = make_positive(rho0)
        rho0 /= np.trace(rho0)
        init_intensity = np.mean(np.multiply(data,1/self.conf["IntensityMap"])) * (rho0.shape[0])
        t0 = density2t(rho0)
        n_t = len(t0)
        # np.multiply(data, self.conf["IntensityMap"])
        # t_to_density(t0)
        t0 = t0 + 0.0001
        t0 = t0 * np.sqrt(init_intensity)

        data = np.real(data)
        data = data.flatten()

        bet = self.conf['Beta']
        # n_data = np.shape(data)[0]


        prediction = np.zeros(m.shape[2]) + 0j

        if bet == 0:
            t = leastsq(self.maxlike_fitness, np.real(t0), args=(data, acc, m, prediction))[0]
            fvalp = np.sum(self.maxlike_fitness(t, data, acc, m, prediction) ** 2)
        else:
            t = \
            leastsq(self.maxlike_fitness_hedged, np.real(t0), args=(data, acc, m, prediction, bet))[0]
            fvalp = np.sum(self.maxlike_fitness_hedged(t, data, acc, m, prediction, bet) ** 2)

        matrix = t_to_density(t)
        intensity = np.trace(matrix)
        rhog = matrix / intensities
        intensity = np.float64(np.real(intensities))

        return [rhog, intensity, fvalp]

    """
    maxlike_fitness(t, data, accidentals, m, prediction)
    Desc: Calculates the diffrence between the current predicted state data and the actual data.

    Parameters
    ----------
    t : ndarray
        T values of the current predicted state.
    data : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The counts of the tomography.
    accidentals : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.   
    m : ndarray with size = (2^numQubits,2^numQubits,number of measurements) 
        The measurements of the tomography in density matrix form.
    prediction : ndarray
        Predicted counts from the predicted state.
    
    Returns
    -------
    val : float
        value of the optimization function.
    """
    def maxlike_fitness(self, t, data, accidentals, m, prediction):

        tm = t_matrix(t)
        # do not change
        rhog = np.dot(tm.conj().transpose(), tm)


        for j in range(len(prediction)):
            prediction[j] = np.float64(np.real(self.conf["IntensityMap"][j] * np.real(np.trace(np.dot(m[:, :, j], rhog))) + accidentals[j]))
            prediction[j] = np.max([prediction[j], 0.01])
        val = (prediction - data) / np.sqrt(prediction)

        val = np.float64(np.real(val))

        return val


    """
    maxlike_fitness_hedged(t, data, accidentals, m, prediction,bet)
    Desc: Calculates the diffrence between the current predicted state data and the actual data using hedged maximum likelihood.

    Parameters
    ----------
    t : ndarray
        T values of the current predicted state.
    data : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The counts of the tomography.
    accidentals : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.   
    m : ndarray with size = (2^numQubits,2^numQubits,number of measurements) 
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
    def maxlike_fitness_hedged(self,t, data, accidentals, m, prediction, bet):

        tm = t_matrix(t)
        rhog = np.dot(tm.conj().transpose(), tm)


        for j in range(len(prediction)):
            prediction[j] = self.conf["IntensityMap"][j] * np.real(np.trace(np.dot(m[:, :, j], rhog))) + accidentals[j]
            prediction[j] = np.max([prediction[j], 0.01])

        hedge = np.repeat(np.real((bet * np.log(np.linalg.det(np.mat(rhog)))) / len(prediction)), len(prediction))
        val = np.sqrt(np.real((((prediction - data) ** 2) / (2 * prediction)) - hedge) + 1000)

        val = np.float64(np.real(val))

        return val


    """
    linear_tomography(data, measurements)
    Desc: Uses linear techniques to find a starting state for maximum likelihood estimation.

    Parameters
    ----------
    data : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The counts of the tomography.
    measurements : ndarray with size = (number of measurements,2*numQubits) 
        The measurements of the tomography in pure state form .
        
    Returns
    -------
    rhog : ndarray with size = (2^numQubits,2^numQubits) 
        The starting predicted state.
    intensity : float
        The predicted overall intensity used to normalize the state.
    """
    def linear_tomography(self,data, measurements, m_set=()):
        if m_set == ():
            m_set = independent_set(measurements)
        if np.isscalar(m_set):
            n = len(data)
            linear_measurements = measurements
            linear_data = data

        else:
            n = np.int(np.sum(m_set))
            linear_measurements = measurements[(np.rot90(m_set == 1.0)[0])]
            linear_data = data[(np.rot90(m_set == 1.0)[0])]

        linear_rhog = np.zeros([measurements.shape[1], measurements.shape[1]])

        b = b_matrix(linear_measurements)
        b_inv = np.linalg.inv(b)

        m = np.zeros([measurements.shape[1], measurements.shape[1], n]) + 0j
        for j in range(n):
            m[:, :, j] = m_matrix(j, linear_measurements, b_inv)
            linear_rhog = linear_rhog + linear_data[j] * m[:, :, j]

        intensities = np.trace(linear_rhog)
        rhog = linear_rhog / intensities

        return [rhog, intensities]

    """
    filter_data(raw_counts, intensities)
    Desc: filters the data into separate arrays

    Parameters
    ----------
    raw_counts : ndarray
        The input data for the current tomography.
    intensities : ndarray with length = number of measurements
        Relative pump power (arb. units) during measurement; used for drift correction.
        
    Returns
    -------
    data : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The counts of the tomography.
    m : ndarray with size = (2^numQubits,2^numQubits,number of measurements) 
        The measurements of the tomography in density matrix form.
    m2 : ndarray with size = (number of measurements,2*numQubits) 
        The measurements of the tomography in pure state form.
    acc : ndarray with length = number of measurements or size =(number of measurements,2^numQubits) for 2 det/qubit
        The singles values of the tomography. Used for accidental correction.
    """
    def filter_data(self, raw_counts, intensities):
        # getting variables
        self.input = raw_counts
        # if not np.isscalar(intensities):
        #     self.conf['DoDriftCorrection'] = 0
        nbits = self.conf['NQubits']
        ndet = self.conf['NDetectors']
        eff = self.conf['Efficiency'][0:2 ** nbits]

        # time values
        t = raw_counts[:, 0]

        # singles values
        n_singles = self.getNumSingles()
        # sings = raw_counts[:, np.arange(n_singles) + 1]
        sings = self.getSingles()

        # state dimension @ qudit size
        # qudit_sizes = self.conf['QuditSizes']
        # if (ndet == 1 and nbits == 1):
        #     qudit_sizes = 2
        # if (ndet == 2):
        #     self.conf['StateDimension'] = 2 ** nbits
        # else:
        #     self.conf['StateDimension'] = np.prod(qudit_sizes)

        # coincidences
        n_coinc = self.getNumCoinc()
        # coinc = raw_counts[:, np.arange(n_singles + 1, n_singles + n_coinc + 1)]
        # if (ndet == 1):
        #     coinc = raw_counts[:, n_singles + 1]
        coinc = self.getCoincidences()

        # settings
        settings = raw_counts[:, np.arange(n_singles + n_coinc + 1, len(raw_counts[0]))]

        # Accidental Correction
        acc = np.zeros_like(coinc)
        if (len(coinc.shape) == 1):
            acc = acc[:, np.newaxis]
        window = self.conf['Window']
        if np.isscalar(window):
            window = np.array([window])
        if any(window) > 0:
            for j in range(n_coinc):
                idx = (multiloop_index(j, 2 * np.ones(nbits)) - 1) * nbits + range(nbits)
                idx = idx.astype(int)
                acc[:, j] = np.prod(np.real(sings[:, idx]), axis=1) * (window[j] * 1e-9 / np.real(t)) ** (nbits - 1)
        if (acc.shape != coinc.shape):
            acc = acc[:, 0]

        # Drift Correction
        self.conf['IntensityMap'] = np.kron(intensities, np.ones(n_coinc))

        # crosstalk
        ctalk = np.array(self.conf['Crosstalk'])[0:2 ** nbits, 0:2 ** nbits]
        crosstalk = ctalk
        if np.ndim(ctalk) >= 3:
            for j in range(ctalk.shape[2]):
                crosstalk[j] = ctalk[:, :, j]

        if (ctalk == []):
            big_crosstalk = np.eye(2 ** nbits)
        else:
            big_crosstalk = crosstalk[:]

        big_crosstalk = big_crosstalk * np.outer(eff, np.ones(n_coinc))

        m = np.zeros([2 ** nbits, 2 ** nbits, np.prod(coinc.shape)]) + 0j
        m2 = np.zeros([np.prod(coinc.shape), 2 ** nbits]) + 0j
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

                m2[j * n_coinc, :] = u[:, 0].conj().transpose()
                for k in range(1):
                    for l in range(2 ** nbits):
                        m[:, :, j + k] = m[:, :, j + k] + m_twiddle[l, :, :] * big_crosstalk[k, l]
                data = coinc

            else:
                for k in range(2 ** nbits):
                    m_twiddle[k, :, :] = np.outer(u[:, k].conj().transpose(), u[:, k])
                    m2[j * n_coinc + k, :] = u[:, k].conj().transpose()
                for k in range(2 ** nbits):
                    for l in range(2 ** nbits):
                        m[:, :, j * (2 ** nbits) + k] = m[:, :, j * (2 ** nbits) + k] + m_twiddle[l, :, :] * \
                                                        big_crosstalk[k, l]

                data = coinc.reshape((np.prod(coinc.shape), 1))
                acc = acc.reshape((np.prod(acc.shape), 1))
        return [data, m, m2, acc]

    ###################
    '''Get Functions'''
    ###################


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
            return self.input[:, np.arange(2*self.conf['NQubits']+1, 2**self.conf['NQubits']+2*self.conf['NQubits']+1)]
        else:
            return self.input[:, self.conf['NQubits']+1]

    """
    getSingles()
    Desc: Returns an array of singles for all the measurments.
    """
    def getSingles(self):
        if (self.conf['NDetectors'] == 2):
            return self.input[:, np.arange(1, 2*self.conf['NQubits']+1)]
        else:
            return self.input[:, np.arange(1, self.conf['NQubits']+1)]

    """
    getTimes()
    Desc: Returns an array of times for all the measurments.
    """
    def getTimes(self):
        return self.tomo_input[:, 0]

    """
        getMeasurements()
        Desc: Returns an array of measurements in pure state form for all the measurments.
        """
    def getMeasurements(self):
        if (self.conf['NDetectors'] == 2):
            return self.input[:, np.arange(2**self.conf['NQubits']+2*self.conf['NQubits']+1, 2**self.conf['NQubits']+4*self.conf['NQubits']+1)]
        else:
            return self.input[:, np.arange(self.conf['NQubits']+2, 3*self.conf['NQubits']+2)]

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
    getBasisMeas(numBits)
    Desc: Returns an array of standard measurments in pure state form for the given number of qubits.
    
    Parameters
    ----------
    numBits : int
        number of qubits you want for each measurement. Default will use the number of qubits in the current configurations.
    """
    def getBasisMeas(self,numBits = -1):
        if(numBits == -1):
            numBits = self.getNumDetPerQubit()
        if(numBits == 1):
            basis = np.array([[1, 0], [0, 1], [(2 ** (-1 / 2)), (2 ** (-1 / 2))], [(2 ** (-1 / 2)), -(2 ** (-1 / 2))],
                              [(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j], [(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j]],
                             dtype=complex)
            m = np.zeros((6 ** self.getNumBits(), 2 * self.getNumBits()), dtype=complex)
            for i in range(m.shape[0]):
                for j in range(0, m.shape[1], 2):
                    bitNumber = np.floor((m.shape[1] - j - 1) / 2)
                    index = int(((i) / 6 ** (bitNumber)) % 6)

                    m[i, j] = basis[index][0]
                    m[i, j + 1] = basis[index][1]
        else:
            basis = np.array([[1, 0],[(2 ** (-1 / 2)), (2 ** (-1 / 2))], [(2 ** (-1 / 2)),(2 ** (-1 / 2)) * 1j]], dtype=complex)
            m = np.zeros((3 ** self.getNumBits(), 2 * self.getNumBits()), dtype=complex)
            for i in range(m.shape[0]):
                for j in range(0, m.shape[1], 2):
                    bitNumber = np.floor((m.shape[1] - j - 1) / 2)
                    index = int(((i) / 3 ** (bitNumber)) % 3)

                    m[i, j] = basis[index][0]
                    m[i, j + 1] = basis[index][1]
        return m


    """
    getTomoInputTemplate()
    Desc: returns a standard template for tomo_input for the given number of qubits.
    
    Parameters
    ----------
    numBits : int
        number of qubits you want for each measurement. Default will use the number of qubits in the current configurations.
    """
    def getTomoInputTemplate(self,numBits = -1):

        measurements = self.getBasisMeas(numBits)

        if(self.getNumDetPerQubit() == 1):
            # For n detectors:
            Tomoinput = np.zeros((6**self.getNumBits(),3*self.getNumBits()+2),dtype=complex)

            # input[:, n_qubit+1]: coincidences
            # coincidences left zeros

            # input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
            Tomoinput[:, np.arange(self.getNumBits() + 2, 3 * self.getNumBits() + 2)] = measurements

        else:
            # For 2n detectors:
            Tomoinput = np.zeros((3**self.getNumBits(),2**self.getNumBits()+4*self.getNumBits()+1),dtype=complex)

            # input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
            # coincidences left zeros

            # input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements
            Tomoinput[:, np.arange(2**self.getNumBits()+2*self.getNumBits()+1, 2**self.getNumBits()+4*self.getNumBits()+1)] = measurements

        # tomo_input[:, 0]: times
        Tomoinput[:,0] = np.ones_like(Tomoinput[:,0])

        # tomo_input[:, np.arange(1, 2 * n_qubit + 1)]: singles
        # singles left zero
        return Tomoinput

    """
    getAllProperties(rho)
    Desc: returns the properties of the given density matrix. 

    Parameters
    ----------
    rho : ndarray with size = (2^numQubits,2^numQubits) 
        The density matrix you would like to know the properties of.
    bounds : boolean
        Set this to true if you want error bounds on your estimated property values. Default is False.
        These are determined with monte carlo simulation.
        
    Returns
    -------
    err_list : ndarray with size = (length of self.err_functions,2)
        The first row is the value of the property and the second is the error bound on the property.
    """
    def getAllProperties(self,rho,bounds = False):
        vals = tomo.getproperties(rho)
        if (bounds > 1):
            err_n = tomo.conf['DoErrorEstimation']
            if(err_n<=1):
                err_n =2
            rhon = tomo.tomography_error_states_generator(err_n)
            [mean, errs, mean_fid, err_fid] = tomo.tomography_error(rho, rhon)
            rho = np.around(rho, 5)
            mean = np.around(mean, 5)
            errs = np.around(errs, 5)
            vals = np.around(vals, 5)
            intensity = round(intensity, 1)
            fval = round(fval, 5)
            err_list = qLib.outputErrorValues(tomo.err_functions, vals, errs, intensity, fval, err_n)
        else:
            outputNames = tomo.err_functions
            outputVals = {}
            for x in range(0, len(outputNames)):
                outputVals[outputNames[x]] = [vals[x],0]
            err_list = [outputNames,outputVals]
        return err_list

    #####################
    '''ERROR FUNCTIONS'''
    #####################

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def tomography_error(self, rho0, rhop):
        if (self.conf['DoErrorEstimation'] == 0):
            self.conf['DoErrorEstimation'] = 2
        n_fun = len(self.err_functions)
        data = np.zeros([self.conf['DoErrorEstimation'], n_fun])
        fid = np.zeros(self.conf['DoErrorEstimation'])
        for j in range(self.conf['DoErrorEstimation']):
            fid[j] = fidelity(rho0, rhop[j, :, :])
            for k in range(n_fun):
                data[j, k] = self.fevel(self.err_functions[k], rhop[j, :, :])

        errors = np.zeros(n_fun)
        means = np.zeros(n_fun)

        for k in range(n_fun):
            errors[k] = np.std(data[:, k])
            means[k] = np.mean(data[:, k])
        mean_f = np.mean(fid)
        err_f = np.std(fid)

        return [means, errors, mean_f, err_f]

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def fevel(self,funcname, *args):
        return eval(funcname)(*args)

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def tomography_error_states_generator(self):
        ndet = self.conf['NDetectors']
        nbits = self.conf['NQubits']
        acc = self.conf['DoAccidentalCorrection']
        rhop = np.zeros([n, 2 ** nbits, 2 ** nbits]) + 0j
        length = len(self.tomo_input[:, 0])
        if ndet == 1:
            time = np.reshape(self.tomo_input[:, 0], (length, 1))
            meas = self.tomo_input[:, np.arange(nbits + 2, len(self.tomo_input[0, :]))]
            for j in range(self.conf['Do']):
                test_data = np.zeros([length, nbits + 1])
                if acc:
                    kk = range(nbits + 1)
                else:
                    kk = np.array([nbits])
                for k in kk:
                    for l in range(length):
                        test_data[l, k] = np.random.poisson(np.real(self.tomo_input[l, k + 1]))

                test_data = np.concatenate((time, test_data, meas), axis=1)

                rhop[j, :, :] = self.state_tomography(test_data, self.intensity)[0]

        elif ndet == 2:
            time = np.reshape(self.tomo_input[:, 0], (length, 1))
            meas = self.tomo_input[:, np.arange(nbits + 2 ** nbits + 1, len(self.tomo_input[0, :]))]
            for j in range(n):
                test_data = np.zeros([length, nbits + 2 ** nbits])
                if acc:
                    kk = range(nbits + 2 ** nbits)
                else:
                    kk = np.arange(nbits, nbits + 2 ** nbits)
                for k in kk:
                    for l in range(length):
                        test_data[l, k] = np.random.poisson(np.real(self.tomo_input[l, k + 1]))

                test_data = np.concatenate((time, test_data, meas), axis=1)

                rhop[j, :, :] = self.state_tomography(test_data, intensities)[0]

        return rhop

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def getproperties(self,rho0):
        return [self.fevel(errf, rho0) for errf in self.err_functions]

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def bellsettings(self,rhog, partsize_init, partsize, t):

        [s, arange, brange, aprange, bprange] = self.bellsettings_range_init(rhog, partsize_init)
        a = 0
        b = 0
        ap = 0
        bp = 0

        for j in range(t):
            [s, a, ap, b, bp, arange, brange, aprange, bprange] = \
                self.bellsettings_range(rhog, partsize, arange, brange, aprange, bprange)

        return [s, a, ap, b, bp]

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def bellsettings_range_init(self,rhog, partsize):
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
                        nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, b + np.pi / 2), rhog)))
                        e_ab = 2 * (npp + nmm) - 1

                        npp = np.real(np.trace(np.dot(coinmat(ap, b), rhog)))
                        nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, b + np.pi / 2), rhog)))
                        e_apb = 2 * (npp + nmm) - 1

                        npp = np.real(np.trace(np.dot(coinmat(a, bp), rhog)))
                        nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, bp + np.pi / 2), rhog)))
                        e_abp = 2 * (npp + nmm) - 1

                        npp = np.real(np.trace(np.dot(coinmat(ap, bp), rhog)))
                        nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, bp + np.pi / 2), rhog)))
                        e_apbp = 2 * (npp + nmm) - 1

                        s = e_ab + e_abp + e_apb + e_apbp - 2 * np.min([e_ab, e_abp, e_apb, e_apbp])

                        if s > sval:
                            sval = s
                            aval = a
                            apval = ap
                            bval = b
                            bpval = bp

        arange_s = [np.max([aval - ((np.pi / 2) / partsize), 0]), np.min([aval + ((np.pi / 2) / partsize), np.pi / 2])]
        aprange_s = [np.max([apval - ((np.pi / 2) / partsize), 0]),
                     np.min([apval + ((np.pi / 2) / partsize), np.pi / 2])]
        brange_s = [np.max([bval - ((np.pi / 2) / partsize), 0]), np.min([bval + ((np.pi / 2) / partsize), np.pi / 2])]
        bprange_s = [np.max([bpval - ((np.pi / 2) / partsize), 0]),
                     np.min([bpval + ((np.pi / 2) / partsize), np.pi / 2])]

        return [sval, arange_s, brange_s, aprange_s, bprange_s]

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def tomography_error_bell(self,rhop, partsize_init, partsize, t, n):
        belldata = np.zeros([n, 5])

        for j in range(n):
            belldata[j, :] = self.bellsettings(rhop[j, :, :], partsize_init, partsize, t)

        bmeans = np.zeros(5)
        berrors = np.zeros(5)

        for m in range(5):
            berrors[m] = np.std(belldata[:, m])
            bmeans[m] = np.mean(belldata[:, m])

        return [berrors, bmeans]

    """
    Function()
    Desc: short desc

    Parameters
    ----------
    x1, x2 : array_like
        Input arrays to be multiplied. If ``x1.shape != x2.shape``, they must be broadcastable to a common shape (which becomes the shape of the output).

    Returns
    -------
    y : ndarray
        The product of `x1` and `x2`, element-wise.
    """
    def bellsettings_range(self, rhog, partsize, arange, brange, aprange, bprange):

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
                        nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, b + np.pi / 2), rhog)))
                        e_ab = 2 * (npp + nmm) - 1

                        npp = np.real(np.trace(np.dot(coinmat(ap, b), rhog)))
                        nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, b + np.pi / 2), rhog)))
                        e_apb = 2 * (npp + nmm) - 1

                        npp = np.real(np.trace(np.dot(coinmat(a, bp), rhog)))
                        nmm = np.real(np.trace(np.dot(coinmat(a + np.pi / 2, bp + np.pi / 2), rhog)))
                        e_abp = 2 * (npp + nmm) - 1

                        npp = np.real(np.trace(np.dot(coinmat(ap, bp), rhog)))
                        nmm = np.real(np.trace(np.dot(coinmat(ap + np.pi / 2, bp + np.pi / 2), rhog)))
                        e_apbp = 2 * (npp + nmm) - 1

                        s = e_ab + e_abp + e_apb + e_apbp - 2 * np.min([e_ab, e_abp, e_apb, e_apbp])

                        if s > sval:
                            sval = s
                            aval = a
                            apval = ap
                            bval = b
                            bpval = bp

        arange_s = [np.max([aval - ((arange[1] - arange[0]) / partsize), 0]),
                    np.min([aval + ((arange[1] - arange[0]) / partsize), np.pi / 2])]
        aprange_s = [np.max([apval - ((aprange[1] - aprange[0]) / partsize), 0]),
                     np.min([apval + ((aprange[1] - aprange[0]) / partsize), np.pi / 2])]
        brange_s = [np.max([bval - ((brange[1] - brange[0]) / partsize), 0]),
                    np.min([bval + ((brange[1] - brange[0]) / partsize), np.pi / 2])]
        bprange_s = [np.max([bpval - ((bprange[1] - bprange[0]) / partsize), 0]),
                     np.min([bpval + ((bprange[1] - bprange[0]) / partsize), np.pi / 2])]

        return [sval, aval, apval, bval, bpval, arange_s, brange_s, aprange_s, bprange_s]
