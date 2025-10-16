from __future__ import print_function

from .TomoFunctions import *
from .TomoDisplay import floatToString
from .TomoClassHelpers import *
from .Utilities import (
    ConfDict,
    NEW_FORMAT_CONFIG_KEY_MAPPING,
    NEW_FORMAT_DATA_KEY_MAPPING,
    OLD_FORMAT_CONFIG_KEYS,
    cast_to_numpy,
    get_raw_measurement_bases_from_data,
    get_all_product_states_from_data,
    get_all_measurements_from_data,
    get_highest_fold_coincidence_count_index,
    parse_np_array,
)
import numpy as np
from scipy.optimize import leastsq
import warnings
from pathlib import Path
import json

from math import modf
import tomllib
import tomli_w
import re

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

"""
    class Tomography()
    Desc: This is the main tomography object that the library is built around. If you are looking to perform multiple tomographies
          it is advised to create multiple Tomography objects.

    See Also
     ------ 
    importEval;importConf;importData
    """


class Tomography:
    # # # # # # # # # # # #
    """Public Variables"""

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
    # conf['ftol']: Relative error desired in the sum of squares. Passed directly to scipy.optimize.leastsq
    # conf['xtol']: Relative error desired in the approximate solution. Passed directly to scipy.optimize.leastsq
    # conf['gtol']: Orthogonality desired between the function vector and the columns of the Jacobian. Passed directly to scipy.optimize.leastsq
    # conf['maxfev']: The maximum number of calls to the function. Passed directly to scipy.optimize.leastsq

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
    """Constructors"""
    # # # # # # # # # #

    """
    Default Constructor
    Desc: This initializes a default tomography object
    """

    def __init__(self, nQ=1):
        self.tomo_input = None
        self.conf = ConfDict(
            [
                ("NQubits", nQ),
                ("NDetectors", 1),
                ("Crosstalk", np.eye(2**nQ)),
                ("Bellstate", 0),
                ("DoDriftCorrection", 0),
                ("DoAccidentalCorrection", 0),
                ("DoErrorEstimation", 0),
                ("Window", [1]),
                ("Efficiency", [1]),
                ("RhoStart", None),
                ("Beta", 0),
                ("ftol", 1.49012e-08),
                ("xtol", 1.49012e-08),
                ("gtol", 0.0),
                ("maxfev", 0),
                ("Method", "MLE"),
            ]
        )
        self.err_functions = ["concurrence", "tangle", "entropy", "linear_entropy", "negativity", "purity"]
        self.mont_carlo_states = list()

    """
    setConfSetting(setting, val)
    Desc: Sets a specific self.conf setting. This is no longer needed but will be kept in the class to
    support old code.

    Parameters
    ----------
    setting : string
        The setting you want to update.
        Possible values are ['NQubits', 'NDetectors', 'Crosstalk', 'Bellstate', 'DoDriftCorrection', 
        'DoAccidentalCorrection', 'DoErrorEstimation', 'Window', 'Efficiency', 'RhoStart', 'Beta','ftol','xtol','gtol'
        'maxfev']
    val : ndarray, int, or string
            The new value you want to the setting to be.
    """

    def setConfSetting(self, setting, val):
        if setting == "get_bell_settings":
            internal_setting = "Bellstate"
        else:
            split = setting.split("_")
            internal_setting = "".join(word.capitalize() for word in split)

        if internal_setting not in self.conf:
            raise ValueError(f"{setting} is not a valid setting for the configuration.")
        self.conf[internal_setting] = val

    """
    importConf(conftxt)
    Desc: Import a text file containing the configuration settings.

    Parameters
    ----------
    conftxt : string
        path to configuration file
    """

    def importConf(self, filename):
        warnings.warn(
            "The old file format (pythoneval.txt, conf.txt, data.txt) will be removed in future versions. This is the last version for which it will be supported."
            "Please switch to using the new configuration (.toml) and data (.json) formats as shown in ExampleFiles. You can import the new files with import_conf()."
        )
        self.import_eval(filename)

    def import_conf(self, filename):
        """Import configuration file with new file format (toml)."""
        with open(Path(filename), "rb") as f:
            conf = tomllib.load(f)
            self._import_conf(conf)

    def _import_conf(self, conf):
        if "get_bell_settings" in conf:
            self.conf["Bellstate"] = conf["get_bell_settings"]

        if "do_drift_correction" in conf:
            self.conf["DoDriftCorrection"] = conf["do_drift_correction"]

        if "do_error_estimation" in conf:
            self.conf["DoErrorEstimation"] = conf["do_error_estimation"]

        if "do_accidental_correction" in conf:
            self.conf["DoAccidentalCorrection"] = conf["do_accidental_correction"]

        if "method" in conf:
            self.conf["Method"] = conf["method"]

        if "starting_matrix" in conf:
            self.conf["RhoStart"] = np.array(conf["starting_matrix"], dtype=complex)

        if "beta" in conf:
            self.conf["beta"] = conf["beta"]

        if "ftol" in conf:
            self.conf["ftol"] = conf["ftol"]

        if "xtol" in conf:
            self.conf["xtol"] = conf["xtol"]

        if "gtol" in conf:
            self.conf["gtol"] = conf["gtol"]

        if "maxfev" in conf:
            self.conf["maxfev"] = conf["maxfev"]

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
        warnings.warn(
            "The old file format (pythoneval.txt, conf.txt, data.txt) will be removed in future versions. This is the last version for which it will be supported."
            "Please switch to using the new configuration (.toml) and data (.json) formats as shown in ExampleFiles. You can import the new files with import_conf()."
        )
        return self.import_eval(datatxt)

    def import_data(self, filename):
        """Import data with new file format (json)."""
        with open(Path(filename), "rb") as f:
            json_dict = json.load(f)
            self._import_data(json_dict)

    def _import_data(self, json_dict):
        for state_name in json_dict["measurement_states"].keys():
            cast_to_numpy(json_dict["measurement_states"], state_name)

        self.conf["measurement_states"] = json_dict["measurement_states"]
        self.conf["data"] = json_dict["data"]

        self.conf["NQubits"] = json_dict["n_qubits"]
        self.conf["NDetectors"] = json_dict["n_detectors_per_qubit"]
        self.conf["NMeasurementsPerQubit"] = json_dict["n_measurements_per_qubit"]
        detector_efficiency_matrix = None
        if "relative_efficiency" in json_dict:
            detector_efficiency_matrix = np.array(json_dict["relative_efficiency"])
        detector_coincidence_window = None
        if "coincidence_window" in json_dict:
            detector_coincidence_window = np.array(json_dict["coincidence_window"])

        # Find which basis states are orthogonal to each other
        orthogonal_bases = {}
        for key1, measurement_1 in json_dict["measurement_states"].items():
            for key2, measurement_2 in json_dict["measurement_states"].items():
                if isStateVector(measurement_1):
                    measurement_1 = np.outer(measurement_1, measurement_1.conj().T)
                    measurement_2 = np.outer(measurement_2, measurement_2.conj().T)
                orthogonal_check = measurement_1 @ measurement_2
                if (orthogonal_check == np.zeros_like(measurement_1)).all() and key1 not in orthogonal_bases:
                    orthogonal_bases[key1] = key2
        singles = []
        coincidences = []
        times = []
        intensities = []

        if self.conf["NDetectors"] > 1:
            measurements = []
            all_meas_projectors = get_all_measurements_from_data(json_dict)
            if "relative_efficiency" in json_dict and not all(
                ["detectors_used" in datum for datum in json_dict["data"]]
            ):
                raise ValueError(
                    "Need to have 'detectors_used' fields filled out in your data array if using more than 1 detector per qubit and you want to account for detector efficiencies!"
                )
            # Find out which measurements occured simultaneously (for normalization)
            simultaneous_measurements = []
            for i, datum1 in enumerate(json_dict["data"]):
                measurement = []
                measurement_basis_1 = [[k, orthogonal_bases[k]] for k in datum1["basis"]]
                measurement_basis_1 = [set(basis) for basis in measurement_basis_1]
                for j, datum2 in enumerate(json_dict["data"]):
                    measurement_basis_2 = [[k, orthogonal_bases[k]] for k in datum2["basis"]]
                    measurement_basis_check = [
                        set(basis) == measurement_basis_1[k] for k, basis in enumerate(measurement_basis_2)
                    ]
                    if all(measurement_basis_check):
                        measurement.append(j)
                if measurement not in simultaneous_measurements:
                    simultaneous_measurements.append(measurement)

            # Construct the measurement basis list, singles and coincidences
            # based on the simultaneous measurements
            for measurement in simultaneous_measurements:
                basis_names = [json_dict["data"][k]["basis"] for k in measurement]
                # Get the projectors corresponding to the measurement from simultaneous_measurements
                measurements.append([all_meas_projectors[k] for k in measurement])
                # Get all the data corresponding to the measurements that happened
                # simultaneously
                measurement_data = [json_dict["data"][k] for k in measurement]
                # Get the index of the highest fold coincidence counts in the data based
                # on the bases defined in the first data point (assuming that all the measurements
                # will have the same structure of counts.
                coincidence_idx = get_highest_fold_coincidence_count_index(len(measurement_data[0]["basis"]))
                # Define expected singles counts array
                singles_counts = np.zeros(self.conf["NDetectors"] * self.conf["NQubits"])
                coincidence_counts = []
                """
                If user defined an intensity, we need to use that to normalize
                (Assuming that the user defined the same relative intensity for
                all the data points corresponding to the same simultaneous measurement.
                This normalizes this simultaneous measurement with the other simultaneous
                measurements, e.g, one simultaneous measurement measures HH, HV, VH, VV
                and has relative intensity 0.9, and another simultaneous measurement measures
                DH, DV, AH, AV and has relative intensity 1. Then, the data points for HH, HV, VH,
                and VV will all have relative_intensity to be defined as 0.9.)
                """
                if "relative_intensity" in measurement_data[0]:
                    measurement_norm = measurement_data[0]["relative_intensity"]
                else:
                    measurement_norm = 1
                # TODO: This will only get the data based on the last simultaneous measurement
                # Is there a case where we can't use only one of the simultaneous measurements
                # to get the efficiency and window?
                efficiencies = []
                window = []
                # Get the singles counts for each detector involved in this measurement
                for datum in measurement_data:
                    # Get the singles for each detector used
                    for count_idx, detector in enumerate(datum["detectors_used"]):
                        singles_counts[detector] = datum["counts"][count_idx]
                    # Get the highest-fold coincidence counts for this measurement
                    coincidence_counts.append(datum["counts"][coincidence_idx])
                    # If detector efficiencies are defined, then we need to use that
                    # for normalization
                    if detector_efficiency_matrix is not None:
                        efficiencies.append(detector_efficiency_matrix[tuple(datum["detectors_used"])])
                    else:
                        efficiencies.append(1)
                    if detector_coincidence_window is not None:
                        window.append(detector_coincidence_window[tuple(datum["detectors_used"])])
                    else:
                        window.append(0)
                if "integration_time" in measurement_data[0]:
                    times.append(float(measurement_data[0]["integration_time"]))
                else:
                    times.append(1.0)
                singles.append(singles_counts)
                coincidences.append(coincidence_counts)

            measurements = np.array(measurements)
            # This measurements array will include elements like (HH, HV, VH, VV)
            # (DH, DV, AH, AV), etc. We only want to keep the first of these, i.e, (HH), (DH), etc.
            # since this is what the current code uses. Essentially, it only wants to know the
            # measurement basis and not the full list of projectors
            self.measurements = measurements[:, 0, :]
            window = np.array(window)
        else:
            efficiencies = np.ones(self.conf["NDetectors"] ** self.conf["NQubits"])
            window = detector_coincidence_window
            for datum in json_dict["data"]:
                idx = get_highest_fold_coincidence_count_index(len(datum["basis"]))
                singles.append(np.array(datum["counts"])[0 : len(datum["basis"])])
                coincidences.append(np.array(datum["counts"])[idx])
                if "integration_time" in datum:
                    times.append(float(datum["integration_time"]))
                else:
                    times.append(1.0)
                if "relative_intensity" in datum:
                    intensities.append(float(datum["relative_intensity"]))
            self.measurements = get_raw_measurement_bases_from_data(json_dict)
        if not len(intensities):
            intensities = np.ones(len(singles))
        elif len(intensities) > 0 and len(intensities) != len(singles):
            raise ValueError("Missing intensities for some measurements!")
        self.time = np.array(times)
        self.singles = np.array(singles, dtype=int)
        self.counts = np.array(coincidences, dtype=int)
        self.intensities = intensities
        self.conf["Efficiency"] = np.array(efficiencies)
        self.conf["Window"] = window
        if "crosstalk" in json_dict:
            input_crosstalk = np.array(json_dict["crosstalk"])
            if input_crosstalk.ndim > 2:
                partial_crosstalk = input_crosstalk[0]
                for i in range(1, input_crosstalk.shape[0]):
                    partial_crosstalk = np.kron(partial_crosstalk, input_crosstalk[i])
                self.conf["Crosstalk"] = partial_crosstalk
            else:
                self.conf["Crosstalk"] = input_crosstalk
        else:
            self.conf["Crosstalk"] = np.eye(2**self.conf["NQubits"])

        # Build tomo input and set it as the last input
        # (in case user doesn't run and only wants to convert config formats)
        self.tomo_input = self.buildTomoInput(
            self.measurements,
            self.counts,
            self.conf["Crosstalk"],
            self.conf["Efficiency"],
            self.time,
            self.singles,
            self.conf["Window"],
            error=self.conf["DoErrorEstimation"],
        )
        self.last_input = np.array(self.tomo_input)
        _ = self.filter_data(self.tomo_input)

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
        warnings.warn(
            "The old file format (pythoneval.txt, conf.txt, data.txt) will be removed in future versions. This is the last version for which it will be supported."
            "Please switch to using the new configuration (.toml) and data (.json) formats as shown in ExampleFiles. You can import the new files with import_conf()."
        )
        return self.import_eval(evaltxt)

    def import_eval(self, filename):
        """Handle old file format in a safe way. importData and importConf both point to here.
        New file format is handled by import_data and import_conf
        """
        tomo_input = None
        intensities = -1
        with open(Path(filename), "r") as f:
            for line in f:
                assignment = line.split("=")
                assignment[1] = assignment[1].strip()
                assignment[0] = assignment[0].strip()

                if assignment[0] == "tomo_input":
                    tomo_input = parse_np_array(assignment[1])
                elif assignment[0] == "intensity":
                    intensities = parse_np_array(assignment[1])
                else:
                    pattern = str(r"""(\[['"][\w]*['"]\])""")
                    match = re.search(pattern, str(assignment[0]))
                    if match:
                        variable_name = match.group(0).strip().strip("[]").strip("'").strip('"').lower()
                        if variable_name == "nqubits":
                            self.conf["NQubits"] = int(assignment[1])

                        elif variable_name == "ndetectors":
                            self.conf["NDetectors"] = int(assignment[1])

                        elif variable_name == "crosstalk":
                            if assignment[1] == "-1":
                                self.conf["Crosstalk"] = np.eye(2**self.conf["NQubits"])
                            else:
                                self.conf["Crosstalk"] = parse_np_array(assignment[1])

                        elif variable_name == "bellstate":
                            if assignment[1] == "yes":
                                self.conf["Bellstate"] = 1
                            else:
                                self.conf["Bellstate"] = 0

                        elif variable_name == "dodriftcorrection":
                            if assignment[1] == "yes":
                                self.conf["DoDriftCorrection"] = 1
                            else:
                                self.conf["DoDriftCorrection"] = 0

                        elif variable_name == "doaccidentalcorrection":
                            if assignment[1] == "yes":
                                self.conf["DoAccidentalCorrection"] = 1
                            else:
                                self.conf["DoAccidentalCorrection"] = 0

                        elif variable_name == "doerrorestimation":
                            self.conf["DoErrorEstimation"] = int(assignment[1])

                        elif variable_name == "efficiency":
                            self.conf["Efficiency"] = parse_np_array(assignment[1])

                        elif variable_name == "rhostart":
                            if assignment[1].strip() == "[]" or assignment[1].strip() == "None":
                                self.conf["RhoStart"] = None
                            else:
                                self.conf["RhoStart"] = parse_np_array(assignment[1])

                        elif variable_name == "usederivative":
                            self.conf["UseDerivative"] = int(assignment[1])

                        elif variable_name == "method":
                            self.conf["Method"] = assignment[1].strip("'").strip('"')

                        elif variable_name == "window":
                            if assignment[1].isdigit():
                                self.conf["Window"] = np.array([int(assignment[1])])
                            else:
                                self.conf["Window"] = np.array(assignment[1])

                        elif variable_name == "beta":
                            self.conf["Beta"] = float(assignment[1])

                        elif variable_name == "ftol":
                            self.conf["ftol"] = float(assignment[1])

                        elif variable_name == "gtol":
                            self.conf["gtol"] = float(assignment[1])

                        elif variable_name == "xtol":
                            self.conf["xtol"] = float(assignment[1])

                        elif variable_name == "maxfev":
                            self.conf["maxfev"] = int(assignment[1])
                        elif variable_name == "nmeasurementsperqubit":
                            self.conf["NMeasurementsPerQubit"] = int(assignment[1])
                        else:
                            raise ValueError(f"{line} is not a valid option for configuration!")

        if tomo_input is not None:
            self.last_input = tomo_input
            self.tomo_input = tomo_input
            self.intensities = intensities
            self.filter_data(tomo_input)
            return self.run_tomography()

    # # # # # # # # # # # # # #
    """Tomography Functions"""
    # # # # # # # # # # # # # #

    """
    StateTomography(measurements, counts)
    Desc: Main function that runs tomography. This function requires a set of measurements and a set of counts.

    Parameters
    ----------
    measurements : ndarray shape = ( Number of measurements , 2*NQubits )
        Each row in the matrix is a set of independent measurements.
    counts : ndarray shape = (Number of measurements, NDetectors**NQubits )
        Each row in the matrix is a set of independent measurements.
    crosstalk : ndarray shape = ( Number of measurements , 2**NQubits,2**NQubits ) (optional)
        The crosstalk matrix to compensate for inefficentcies in your beam splitter.
    efficiency : 1darray with length = NQubits*NDetectors (optional)
        The relative efficienies between your detector pairs.
    time :  1darray with length = Number of measurements (optional)
        The total duration of each measurment.
    singles : ndarray shape = ( Number of measurements , NDetectors*NQubits ) (optional)
        The singles counts on each detector.
    window : 1darray with length = NQubits*NDetectors (optional)
        The coincidence window duration for each detector pair.
    error : int (optional)
        The number of monte carlo states used to estimate the properties of the state.
    intensities : 1darray with length = Number of measurements (optional)
        The relative intensity of each measurement. Used for drift correction.
    method : ['MLE','HMLE','BME','LINER'] (optional)
        Which method to use to run tomography. Default is MLE

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

    def StateTomography(
        self,
        measurements,
        counts,
        crosstalk=-1,
        efficiency=0,
        time=-1,
        singles=-1,
        window=0,
        error=0,
        intensities=-1,
        method="MLE",
    ):
        tomo_input = self.buildTomoInput(measurements, counts, crosstalk, efficiency, time, singles, window, error)
        return self.StateTomography_Matrix(tomo_input, intensities, method=method)

    def run_tomography(self):
        tomo_input = self.tomo_input

        if tomo_input is None:
            tomo_input = self.buildTomoInput(
                self.measurements,
                self.counts,
                self.conf["Crosstalk"],
                self.conf["Efficiency"],
                self.time,
                self.singles,
                self.conf["Window"],
                error=self.conf["DoErrorEstimation"],
            )
        return self.StateTomography_Matrix(tomo_input, self.intensities, method=self.conf["Method"])

    """
    StateTomography_Matrix(tomo_input, intensities)
    Desc: Main function that runs tomography.

    Parameters
    ----------
    tomo_input : ndarray
        The input data for the current tomography.
    intensities : 1darray with length = number of measurements (optional)
        Relative pump power (arb. units) during measurement; used for drift correction. Default will be an array of ones
    method : ['MLE','HMLE','BME','LINER'] (optional)
            Which method to use to run tomography. Default is MLE

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

    def StateTomography_Matrix(self, tomo_input, intensities=-1, method="MLE", _saveState=True):
        # define a uniform intensity if not stated
        if intensities is None or isinstance(intensities, int):
            self.intensities = np.ones(tomo_input.shape[0])
        elif len(intensities.shape) == 1 and intensities.shape[0] == tomo_input.shape[0]:
            self.intensities = intensities
        else:
            raise ValueError("Invalid intensities array")
        if any(self.intensities - np.ones_like(self.intensities)) > 10**-6:
            self.conf["DoDriftCorrection"] = 1

        # filter the data
        self.last_input = tomo_input
        self.conf["Method"] = method
        [coincidences, measurements_densities, measurements_pures, accidentals, overall_norms] = self.filter_data(
            tomo_input
        )

        # get the starting state from tomography_LINEAR if not defined
        starting_matrix = self.conf["RhoStart"]
        if method.upper() == "LINEAR" or not isinstance(starting_matrix, np.ndarray):
            try:
                [starting_matrix, inten_linear] = self.tomography_LINEAR(
                    coincidences, measurements_pures, overall_norms
                )
                # Currently linear tomography gets the phase wrong. So a temporary fix is to just transpose it.
                starting_matrix = starting_matrix.transpose()
                starting_matrix = make_positive(starting_matrix)
                starting_matrix = starting_matrix / np.trace(starting_matrix)
            except:
                raise RuntimeError("Failed to run linear Tomography")

        # Run tomography and find an estimate for the state
        if method == "MLE":
            [rhog, intensity, fvalp] = self.tomography_MLE(
                starting_matrix, coincidences, measurements_densities, accidentals, overall_norms
            )
        elif method.upper() == "HMLE":
            [rhog, intensity, fvalp] = self.tomography_HMLE(
                starting_matrix, coincidences, measurements_densities, accidentals, overall_norms
            )
        # elif method.upper() == "BME":
        #     [rhog, intensity, fvalp] = self.tomography_BME(starting_matrix, coincidences, measurements_densities,accidentals,overall_norms)
        elif method.upper() == "LINEAR":
            [rhog, intensity, fvalp] = [starting_matrix, inten_linear, 0]
        else:
            raise ValueError("Invalid Method name: " + str(method))

        if _saveState:
            # save the results
            self.last_rho = rhog.copy()
            self.last_intensity = intensity
            self.last_fval = fvalp
            # Reset the monte carlo states after doing a tomography
            self.mont_carlo_states = list([[rhog, intensity, fvalp]])

        return [rhog, intensity, fvalp]

    # This is a temporary function. It's here in case there is old code that still uses state_tomography.
    # Eventually this should be removed in a later version
    def state_tomography(self, tomo_input, intensities=-1, method="MLE", _saveState=True):
        warnings.warn(
            "state_tomography will be removed in a future version. It is replaced by StateTomography_Matrix which does the same exact thing.",
            DeprecationWarning,
        )
        return self.StateTomography_Matrix(tomo_input, intensities, method, _saveState)

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

    def tomography_MLE(self, starting_matrix, coincidences, measurements, accidentals, overall_norms=-1):
        # If overall_norms not given then assume uniform
        if not isinstance(overall_norms, np.ndarray):
            overall_norms = np.ones(coincidences.shape[0])
        elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
            raise ValueError("Invalid intensities array")

        init_intensity = np.mean(np.multiply(coincidences, 1 / overall_norms)) * (starting_matrix.shape[0])
        starting_tvals = density2t(starting_matrix)
        starting_tvals = starting_tvals + 0.0001
        starting_tvals = starting_tvals * np.sqrt(init_intensity)

        coincidences = np.real(coincidences)
        coincidences = coincidences.flatten()

        final_tvals = leastsq(
            maxlike_fitness,
            np.real(starting_tvals),
            args=(coincidences, accidentals, measurements, overall_norms),
            ftol=self.conf["ftol"],
            xtol=self.conf["xtol"],
            gtol=self.conf["gtol"],
            maxfev=self.conf["maxfev"],
        )[0]
        fvalp = np.sum(maxlike_fitness(final_tvals, coincidences, accidentals, measurements, overall_norms) ** 2)

        final_matrix = t_to_density(final_tvals, normalize=False)
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

    def tomography_HMLE(self, starting_matrix, coincidences, measurements, accidentals, overall_norms=-1):
        # If overall_norms not given then assume uniform
        if not isinstance(overall_norms, np.ndarray):
            overall_norms = np.ones(coincidences.shape[0])
        elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
            raise ValueError("Invalid intensities array")

        init_intensity = np.mean(np.multiply(coincidences, 1 / overall_norms)) * (starting_matrix.shape[0])
        starting_tvals = density2t(starting_matrix)
        starting_tvals = starting_tvals + 0.0001
        starting_tvals = starting_tvals * np.sqrt(init_intensity)

        coincidences = np.real(coincidences)
        coincidences = coincidences.flatten()

        bet = self.conf["Beta"]
        if bet > 0:
            final_tvals = leastsq(
                maxlike_fitness_hedged,
                np.real(starting_tvals),
                args=(coincidences, accidentals, measurements, bet, overall_norms),
                ftol=self.conf["ftol"],
                xtol=self.conf["xtol"],
                gtol=self.conf["gtol"],
                maxfev=self.conf["maxfev"],
            )[0]
            fvalp = np.sum(
                maxlike_fitness_hedged(final_tvals, coincidences, accidentals, measurements, bet, overall_norms) ** 2
            )
        else:
            raise ValueError("To use Hedged Maximum Likelihood, Beta must be a positive number.")

        final_matrix = t_to_density(final_tvals, normalize=False)
        intensity = np.trace(final_matrix)
        final_matrix = final_matrix / np.trace(final_matrix)

        intensity = np.float64(np.real(intensity))

        return [final_matrix, intensity, fvalp]

    # todo: change these back to double quotes when it is ready to be displayed on documentation page.
    """
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
    """
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

    def tomography_LINEAR(self, coincidences, measurements, overall_norms=-1):
        # If overall_norms not given then assume uniform
        if not isinstance(overall_norms, np.ndarray):
            overall_norms = np.ones(coincidences.shape[0])
        elif not (len(overall_norms.shape) == 1 and overall_norms.shape[0] == coincidences.shape[0]):
            raise ValueError("Invalid intensities array")

        coincidences = coincidences.flatten()

        pauli_basis = generalized_pauli_basis(self.getNumQubits())
        stokes_measurements = (
            np.array([get_stokes_parameters(m, pauli_basis) for m in measurements]) / 2 ** self.getNumQubits()
        )
        freq_array = coincidences / overall_norms

        B_inv = np.linalg.inv(np.matmul(stokes_measurements.T, stokes_measurements))
        stokes_params = np.matmul(stokes_measurements.T, freq_array)
        stokes_params = np.matmul(B_inv, stokes_params)
        linear_rhog = np.multiply(pauli_basis, stokes_params[:, np.newaxis, np.newaxis])
        linear_rhog = np.sum(linear_rhog, axis=0)

        intensity = np.trace(linear_rhog)
        rhog = linear_rhog / intensity

        return [rhog, intensity]

    """
    filter_data(tomo_input)
    Desc: Filters the data into separate arrays.

    Parameters
    ----------
    tomo_input : ndarray
        The input data for the current tomography.

    Returns
    -------
    coincidences : ndarray shape = (Number of measurements, NDetectors**NQubits)
        The counts of the tomography.
    measurements_densities : ndarray with shape = (NDetectors*number of measurements,2^numQubits, 2^numQubits) 
        The measurements of the tomography in density matrix form.
    measurements_pures : ndarray with shape = (NDetectors*number of measurements, 2^numQubits) 
        The measurements of the tomography in pure state form.
    acc : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The accidental values of the tomography. Used for accidental correction.
    overall_norms : 1darray with length = number of measurements or length = number of measurements * 2^numQubits for 2 det/qubit
        The relative weights of each measurment. Used for drift correction.
    """

    def filter_data(self, tomo_input):
        # Check this settings
        self.checkForInvalidSettings()

        # getting variables
        nbits = self.conf["NQubits"]
        times = self.getTimes()
        sings = self.getSingles()
        coinc = self.getCoincidences()
        n_coinc = self.getNumCoinc()
        window = self.conf["Window"]

        eff = self.conf["Efficiency"]
        crosstalk = self.conf["Crosstalk"]
        if isinstance(eff, int):
            eff = np.ones(self.conf["NDetectors"] ** self.conf["NQubits"])
            self.conf["Efficiency"] = eff

        overall_norms = np.kron(self.intensities, eff)

        # Accidental Correction
        acc = np.zeros_like(coinc)
        if self.conf["DoAccidentalCorrection"] == 1:
            scalerIndex = np.concatenate((np.ones(nbits - 2), [2, 2]))
            additiveIndex = np.array([0, 1])
            for j in range(2, nbits):
                additiveIndex = np.concatenate(([2 * j], additiveIndex))
            if len(coinc.shape) == 1:
                acc = acc[:, np.newaxis]
            for j in range(n_coinc):
                index = bin(j).split("b")[1]
                index = "0" * (nbits - len(index)) + index
                index = [int(char) for char in index]
                index = index * scalerIndex + additiveIndex
                index = np.array(index, dtype=int)
                acc[:, j] = np.prod(np.real(sings[:, tuple(index)]), axis=1) * (window[j] * 1e-9 / np.real(times)) ** (
                    nbits - 1
                )
            if acc.shape != coinc.shape:
                acc = acc[:, 0]

        # Get measurements
        measurements_raw = self.getMeasurements()
        measurements_pures = np.zeros([np.prod(coinc.shape), 2**nbits], dtype=complex)
        # NEW WAY
        measurements_densities = np.zeros((tomo_input.shape[0], self.getNumCoinc(), 2**nbits, 2**nbits), dtype=complex)
        for j in range(tomo_input.shape[0]):
            meas_basis_densities = np.zeros([2**nbits, 2**nbits, 2**nbits]) + 0j
            meas_basis_pures = 1
            for k in range(nbits):
                alpha = measurements_raw[j][2 * k]
                beta = measurements_raw[j][2 * k + 1]
                psi_transmit = np.array([alpha, beta])
                psi_reflect = np.array([np.conj(beta), np.conj(-alpha)])
                meas_pure = np.outer((np.array([1, 0])), psi_transmit) + np.outer((np.array([0, 1])), psi_reflect)
                # Changed from tensor_product to np.kron
                meas_basis_pures = np.kron(meas_basis_pures, meas_pure)
            for k in range(2**nbits):
                meas_basis_densities[k, :, :] = np.outer(
                    meas_basis_pures[:, k].conj().transpose(), meas_basis_pures[:, k]
                )
            for k in range(self.getNumCoinc()):
                for l in range(2**nbits):
                    measurements_pures[j * n_coinc + k, :] = meas_basis_pures[:, k].conj().transpose()
                    measurements_densities[j, k, :, :] = (
                        measurements_densities[j, k, :, :] + meas_basis_densities[l, :, :] * crosstalk[k, l]
                    )
        # Flatten the arrays
        measurements_densities = measurements_densities.reshape(
            (
                np.prod(measurements_densities.shape[:2]),
                measurements_densities.shape[2],
                measurements_densities.shape[3],
            )
        )

        # Normalize measurements
        for j in range(measurements_pures.shape[0]):
            norm = np.dot(measurements_pures[j].conj(), measurements_pures[j])
            measurements_pures[j] = measurements_pures[j] / np.sqrt(norm)

        for j in range(measurements_densities.shape[0]):
            norm = np.trace(measurements_densities[j])
            measurements_densities[j] = measurements_densities[j] / norm

        # Flatten the arrays
        coincidences = coinc.reshape((np.prod(coinc.shape), 1))
        acc = acc.reshape((np.prod(acc.shape), 1))
        return [coincidences, measurements_densities, measurements_pures, acc, overall_norms]

    """
    buildTomoInput(tomo_input, intensities)
    desc : This function build an input matrix based on a variety of inputs.

    Parameters
    ----------
    measurements : ndarray shape = ( Number of measurements , 2*NQubits )
        Each row in the matrix is a set of independent measurements.
    counts : ndarray shape = (Number of measurements, NDetectors**NQubits )
        Each row in the matrix is a set of independent measurements.
    crosstalk : ndarray shape = ( Number of measurements , 2**NQubits,2**NQubits ) (optional)
        The crosstalk matrix to compensate for inefficentcies in your beam splitter.
    efficiency : 1darray with length = NQubits*NDetectors (optional)
        The relative efficienies between your detector pairs.
    time :  1darray with length = Number of measurements (optional)
        The total duration of each measurment.
    singles : ndarray shape = ( Number of measurements , 2*NQubits ) (optional)
        The singles counts on each detector.
    window : 1darray with length = NQubits*NDetectors (optional)
        The coincidence window duration for each detector pair.
    error : int (optional)
        The number of monte carlo states used to estimate the properties of the state.
    intensities : 1darray with length = Number of measurements (optional)
        The relative intensity of each measurement. Used for drift correction.
    method : ['MLE','HMLE','BME','LINER'] (optional)
        Which method to use to run tomography. Default is MLE

    Returns
    -------
    tomo_input : ndarray
        The input data for the current tomography.
    intensities : 1darray with length = number of measurements
        Relative pump power (arb. units) during measurement; used for drift correction. Default will be an array of ones
    """

    def buildTomoInput(self, measurements, counts, crosstalk, efficiency, time, singles, window, error):
        ################
        # measurements #
        ################
        # Get the number of qubits based on the measurement matrix
        if len(measurements.shape) != 2 or measurements.shape[1] % 2 != 0:
            raise ValueError("Invalid measurements matrix")
        else:
            self.conf["NQubits"] = int(measurements.shape[1] / 2)
        # TODO: make sure measurements span the entire space

        ##########
        # counts #
        ##########
        # Check if counts has right dimensions
        if counts.shape[0] != measurements.shape[0]:
            raise ValueError("Number of Counts does not match the number of measurements")

        # Determine if 2det/qubit from counts matrix
        try:
            if counts.shape[1] == 1:
                self.conf["NDetectors"] = 1
            elif counts.shape[1] == self.conf["NQubits"] * 2:
                self.conf["NDetectors"] = 2
            else:
                raise ValueError(
                    "The second axis of counts does not have the right dimension. Should be 1 or 2*NQubits for 2det"
                )
        except:
            self.conf["NDetectors"] = 1

        ##############
        # efficiency #
        ##############
        if (
            not isinstance(efficiency, int)
            and self.conf["NDetectors"] > 1
            and efficiency.shape != (self.conf["NDetectors"] ** self.conf["NQubits"],)
        ):
            raise ValueError("Invalid efficiency array. Needs to be an array of size (n_detectors**n_qubits,)")
        else:
            self.conf["Efficiency"] = efficiency

        ##########
        # window #
        ##########
        if (
            not isinstance(window, int)
            and self.conf["NDetectors"] > 1
            and window.shape != (self.conf["NDetectors"] ** self.conf["NQubits"],)
        ):
            raise ValueError("Invalid window array. Needs to be an array of size (n_detectors**n_qubits,)")
        else:
            self.conf["Window"] = window

        ################
        # time/singles #
        ################
        if (not isinstance(time, int)) or (not isinstance(singles, int)) or (not isinstance(window, int)):
            # Check if time has right length
            if isinstance(time, int):
                time = np.ones(measurements.shape[0])
            elif not (len(time.shape) == 1 and time.shape[0] == measurements.shape[0]):
                ValueError("Invalid time array")
            # Check if singles has right dimensions
            if isinstance(singles, int):
                if self.conf["NDetectors"] == 1:
                    singles = np.zeros((measurements.shape[0], self.conf["NQubits"]))
                elif self.conf["NDetectors"] == 2:
                    singles = np.zeros((measurements.shape[0], 2 * self.conf["NQubits"]))
            else:
                if self.conf["NDetectors"] == 1 and self.conf["NQubits"] == 1:
                    if singles.shape != (measurements.shape[0],) and singles.shape != (measurements.shape[0], 1):
                        raise ValueError("Invalid singles matrix")
                    # if the singles vector has the form (x,), this changes it to (x,1)
                    if np.ndim(singles) == 1:
                        singles = np.atleast_2d(singles).transpose()
                elif self.conf["NDetectors"] == 1 and self.conf["NQubits"] != 1:
                    if singles.shape != (measurements.shape[0], self.conf["NQubits"]):
                        raise ValueError("Invalid singles matrix")
                elif self.conf["NDetectors"] == 2:
                    if singles.shape != (measurements.shape[0], 2 * self.conf["NQubits"]):
                        raise ValueError("Invalid Singles matrix")
            # Check if window has right length
            if isinstance(window, int):
                self.conf["Window"] = window
            elif len(window.shape) == 1 and window.shape[0] == self.getNumCoinc():
                self.conf["Window"] = window
            else:
                raise ValueError("Invalid window array")
        else:
            time = np.ones(measurements.shape[0])
            singles = np.zeros((measurements.shape[0], self.getNumSingles()))
            self.conf["Window"] = window

        #############
        # crosstalk #
        #############
        if not isinstance(crosstalk, int) and crosstalk.shape != (
            2 ** self.conf["NQubits"],
            2 ** self.conf["NQubits"],
        ):
            raise ValueError("Invalid crosstalk matrix")
        else:
            self.conf["Crosstalk"] = crosstalk

        #########
        # error #
        #########
        self.conf["DoErrorEstimation"] = error

        ##############
        # tomo_input #
        ##############
        # Here we build the tomo_input matrix and then return this to be fed to StateTomography_Matrix
        n_qubit = self.conf["NQubits"]
        if self.conf["NDetectors"] == 1:
            tomo_input = np.zeros((measurements.shape[0], 3 * n_qubit + 2), dtype=complex)
            # times
            tomo_input[:, 0] = time
            # singles
            tomo_input[:, 1 : n_qubit + 1] = singles
            # counts
            tomo_input[:, n_qubit + 1] = counts
            # measurements
            tomo_input[:, n_qubit + 2 : 3 * n_qubit + 2] = measurements
        else:
            tomo_input = np.zeros((measurements.shape[0], 2**n_qubit + 4 * n_qubit + 1), dtype=complex)
            # times
            tomo_input[:, 0] = time
            # singles
            tomo_input[:, 1 : 2 * n_qubit + 1] = singles
            # counts
            tomo_input[:, 2 * n_qubit + 1 : 2**n_qubit + 2 * n_qubit + 1] = counts
            # measurements
            tomo_input[:, 2**n_qubit + 2 * n_qubit + 1 : 2**n_qubit + 4 * n_qubit + 1] = measurements

        return tomo_input

    """
    checkForInvalidSettings()
    desc : Checks the current settings and throws errors if any are invalid.
    """

    def checkForInvalidSettings(self):
        # Check Accidentals and window
        if self.conf["nqubits"] == 1:
            if self.conf["DoAccidentalCorrection"]:
                raise ValueError(
                    "Invalid Conf settings. Accidental Correction can not be done for single qubit tomography"
                )
        if self.conf["DoAccidentalCorrection"]:
            try:
                win = np.array(self.conf["window"], dtype=float)
                if len(win.shape) != 1 or len(win) != self.getNumCoinc():
                    raise
                self.conf["window"] = win
            except:
                raise ValueError(
                    "Invalid Conf settings. Window should have length "
                    + str(self.getNumCoinc())
                    + " with the given settings."
                )
        # Efficicieny and Ndetectors
        if not self.conf["NDetectors"] in [1, 2]:
            raise ValueError(
                "Invalid Conf settings. NDetectors can be either 1 or 2, corresponding to the number of detectors per qubit."
            )
        elif self.conf["NDetectors"] == 2:
            try:
                eff = np.array(self.conf["Efficiency"], dtype=float)
                if isinstance(self.conf["Efficiency"], int):
                    eff = np.ones(self.conf["NDetectors"] ** self.conf["NQubits"])
                if len(eff.shape) != 1 or len(eff) != self.getNumCoinc():
                    raise
            except:
                raise ValueError(
                    "Invalid Conf settings. Efficiency should have length "
                    + str(self.getNumCoinc())
                    + " with the given settings."
                )
        elif self.conf["NDetectors"] == 1:
            eff = np.ones(self.conf["NDetectors"] ** self.conf["NQubits"])
            self.conf["Efficiency"] = eff
        # Crosstalk
        try:
            correctSize = int(np.floor(2 ** self.getNumQubits() + 0.01))
            c = self.conf["crosstalk"]
            if isinstance(c, int):
                c = np.eye(correctSize)
            else:
                c = np.array(c)

            # make sure it has the right shape
            if len(c.shape) != 2 or c.shape[0] != c.shape[1]:
                raise ValueError(
                    "Invalid Conf settings. Crosstalk should be an array with shape ("
                    + str(correctSize)
                    + ","
                    + str(correctSize)
                    + ") with the given settings."
                )
            # make sure it has the right size
            if c.shape[0] != correctSize:
                if np.any(c - np.eye(c.shape[0]) > 10**-6):
                    # Throw error if not the right size
                    raise ValueError(
                        "Invalid Conf settings. Crosstalk should be an array with shape ("
                        + str(correctSize)
                        + ","
                        + str(correctSize)
                        + ") with the given settings."
                    )
                else:
                    # If the given cross talk is an identity matrix just fix the size.
                    c = np.eye(correctSize)
            self.conf["crosstalk"] = c
        except:
            raise ValueError(
                "Invalid Conf settings. Crosstalk should be an array with shape ("
                + str(correctSize)
                + ","
                + str(correctSize)
                + ") with the given settings."
            )

        # Check starting matrix
        starting_matrix = self.conf["RhoStart"]
        if starting_matrix is not None and np.array(starting_matrix).shape != (
            2 ** self.conf["NQubits"],
            2 ** self.conf["NQubits"],
        ):
            raise ValueError("Starting matrix needs to be of shape (2**n_qubits, 2**n_qubits)")

    # # # # # # # # # #
    """Get Functions"""
    # # # # # # # # # #

    """
    getCoincidences()
    Desc: Returns an array of counts for all the measurements.
    """

    def getCoincidences(self):
        if self.conf["NDetectors"] == 2:
            return np.real(
                self.last_input[
                    :, np.arange(2 * self.conf["NQubits"] + 1, 2 ** self.conf["NQubits"] + 2 * self.conf["NQubits"] + 1)
                ]
            )
        else:
            return np.real(self.last_input[:, self.conf["NQubits"] + 1])

    """
    getSingles()
    Desc: Returns an array of singles for all the measurements.
    """

    def getSingles(self):
        if self.conf["NDetectors"] == 2:
            return np.real(self.last_input[:, np.arange(1, 2 * self.conf["NQubits"] + 1)])
        else:
            return np.real(self.last_input[:, np.arange(1, self.conf["NQubits"] + 1)])

    """
    getTimes()
    Desc: Returns an array of times for all the measurements.
    """

    def getTimes(self):
        return self.last_input[:, 0]

    """
    getMeasurements()
    Desc: Returns an array of measurements in pure state form for all the measurements.
    """

    def getMeasurements(self):
        if self.conf["NDetectors"] == 2:
            return self.last_input[
                :,
                np.arange(
                    2 ** self.conf["NQubits"] + 2 * self.conf["NQubits"] + 1,
                    2 ** self.conf["NQubits"] + 4 * self.conf["NQubits"] + 1,
                ),
            ]
        else:
            return self.last_input[:, np.arange(self.conf["NQubits"] + 2, 3 * self.conf["NQubits"] + 2)]

    """
    getNumQubits()
    Desc: returns the number of qubits for the current configurations.
    """

    def getNumQubits(self):
        return self.conf["NQubits"]

    """
    getNumCoinc()
    Desc: Returns the number of coincidences per measurement for the current configurations.
    """

    def getNumCoinc(self):
        if self.conf["NDetectors"] == 2:
            return 2 ** self.conf["NQubits"]
        else:
            return 1

    """
    getNumSingles()
    Desc: Returns the number of singles per measurement for the current configurations.
    """

    def getNumSingles(self):
        return self.conf["NDetectors"] * self.conf["NQubits"]

    """
    getNumDetPerQubit()
    Desc: returns the number of detectors per qubit for the current configurations.
    """

    def getNumDetPerQubit(self):
        return self.conf["NDetectors"]

    """
    getStandardBasis(numBit,numDet = -1)
    Desc: Returns an array of standard measurements in separated pure state form for the given number of qubits. It is the same format used in the tomo_input matrix.

    Parameters
    ----------
    numBits : int
        number of qubits you want for each measurement. Default will use the number of qubits in the current configurations.
    numDet : 1 or 2
        Number of detectors for each measurement Default will use the number of qubits in the current configurations.   
    """

    def getStandardBasis(self, numBits=-1, numDet=-1):
        # check if numBits is an int
        if not (isinstance(numBits, int)):
            raise ValueError("numBits must be an integer")
        if numBits < 1:
            numBits = self.getNumQubits()

        # check if numDet is not 1 or 2
        if not (numDet in [-1, 1, 2]):
            raise ValueError("numDet must be an 1 or 2")
        if numDet < 1:
            numDet = self.getNumDetPerQubit()

        # Define start meas basis
        if numDet == 1:
            basis = np.array(
                [
                    [1, 0],
                    [0, 1],
                    [(2 ** (-1 / 2)), (2 ** (-1 / 2))],
                    [(2 ** (-1 / 2)), -(2 ** (-1 / 2))],
                    [(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j],
                    [(2 ** (-1 / 2)), -(2 ** (-1 / 2)) * 1j],
                ],
                dtype=complex,
            )
        else:
            basis = np.array(
                [[1, 0], [(2 ** (-1 / 2)), (2 ** (-1 / 2))], [(2 ** (-1 / 2)), (2 ** (-1 / 2)) * 1j]], dtype=complex
            )
        numBasis = basis.shape[0]
        m = np.zeros((numBasis**numBits, 2 * numBits), dtype=complex)
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

    Parameters
    ----------
    numBits : int
        number of qubits you want for each measurement. Default will use the number of qubits in the current configurations.
    numDet : 1 or 2
        Number of detectors for each measurement Default will use the number of qubits in the current configurations.   

    Returns
    -------
    Tomoinput : ndarray
        The input data for the current tomography.
    """

    def getTomoInputTemplate(self, numBits=-1, numDet=-1):
        # check if numBits is an int
        if not (isinstance(numBits, int)):
            raise ValueError("numBits must be an integer")
        if numBits < 1:
            numBits = self.getNumQubits()

        # check if numDet is not 1 or 2
        if not (numDet in [-1, 1, 2]):
            raise ValueError("numDet must be an 1 or 2")
        if numDet < 1:
            numDet = self.getNumDetPerQubit()

        measurements = self.getStandardBasis(numBits, numDet)

        if numDet == 1:
            # For n detectors:
            Tomoinput = np.zeros((6**numBits, 3 * numBits + 2), dtype=complex)

            # input[:, n_qubit+1]: coincidences
            # coincidences left zeros

            # input[:, np.arange(n_qubit+2, 3*n_qubit+2)]: measurements
            Tomoinput[:, np.arange(numBits + 2, 3 * numBits + 2)] = measurements

        else:
            # For 2n detectors:
            Tomoinput = np.zeros((3**numBits, 2**numBits + 4 * numBits + 1), dtype=complex)

            # input[:, np.arange(2*n_qubit+1, 2**n_qubit+2*n_qubit+1)]: coincidences
            # coincidences left zeros

            # input[:, np.arange(2**n_qubit+2*n_qubit+1, 2**n_qubit+4*n_qubit+1)]: measurements
            Tomoinput[:, np.arange(2**numBits + 2 * numBits + 1, 2**numBits + 4 * numBits + 1)] = measurements

        # tomo_input[:, 0]: times
        Tomoinput[:, 0] = np.ones_like(Tomoinput[:, 0])

        # tomo_input[:, np.arange(1, 2 * n_qubit + 1)]: singles
        # singles left zero
        return Tomoinput

    # # # # # # # # # # #
    """ERROR FUNCTIONS"""
    # # # # # # # # # # #

    """
    getProperties(bounds)
    Desc: Returns all the properties of the given density matrix.

    Parameters
    ----------
    bounds : int (optional)
        The number of monte carlo runs you want to perform to get a better estimate of each property. Default will use whatever is set in the conf settings.

    Returns
    -------
    vals : ndarray with shape = (length of self.err_functions, 2)
        The first col is the name of the property.
        The second col is the value of the property.
        The third col is the error bound on the property.
    """

    def getProperties(self, bounds=-1):
        # if bounds not set use the conf settings
        if bounds == -1:
            bounds = self.conf["DoErrorEstimation"]

        # generate states if needed
        if bounds > len(self.mont_carlo_states) - 1:
            self.tomography_states_generator(bounds - len(self.mont_carlo_states) + 1)
        states = np.array(self.mont_carlo_states, dtype="O")[: bounds + 1, :]
        if states.shape[0] == 1:
            vals1 = np.array(
                [["intensity", np.mean(states[:, 1]), "NA"], ["fval", np.mean(states[:, 2]), "NA"]], dtype="O"
            )
        else:
            vals1 = np.array(
                [
                    ["intensity", np.mean(states[:, 1]), np.std(states[:, 1], ddof=1)],
                    ["fval", np.mean(states[:, 2]), np.std(states[:, 2], ddof=1)],
                ],
                dtype="O",
            )
        vals2 = getProperties_helper_bounds(self.err_functions, states[:, 0])

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
    intenp : 1darray with length = n
        The intensity of each approximate density matrices.
    fvalp : 1darray with length = n
        The fval of the regression associated with each approximate density matrices.
    """

    def tomography_states_generator(self, n):
        # Save the last data so we can restore it later
        last_input = self.last_input.copy()

        ndet = self.conf["NDetectors"]
        acc = self.conf["DoAccidentalCorrection"]
        time = self.getTimes()
        meas = self.getMeasurements()
        singles = self.getSingles()
        counts = self.getCoincidences()

        for j in range(n):
            # Re-sample counts and singles
            test_counts = np.random.poisson(counts)
            test_singles = np.random.poisson(singles)
            if len(test_counts.shape) == 1:
                test_counts = np.array([test_counts]).T
            test_data = np.concatenate((np.array([time]).T, test_singles, test_counts, meas), axis=1)
            [rhop, intenp, fvalp] = self.StateTomography_Matrix(
                test_data, self.intensities, method=self.conf["method"], _saveState=False
            )
            self.mont_carlo_states.append([rhop, intenp, fvalp])
        # Restore the last tomo_input matrix to the original one
        self.last_input = last_input
        return self.mont_carlo_states

    """
    getBellSettings(bounds)
    Desc: Returns the optimal measurment settings for the CHSH bell inequality. In Progress, have not checked.
    """

    def getBellSettings(self, partsize_init=9, partsize=5, t=3, bounds=-1):
        if self.conf["nqubits"] != 2:
            raise ValueError("Invalid Conf settings. Bell state can only be found for 2 qubit tomography")
        rho = self.last_rho
        # if bounds not set use the conf settings
        if bounds == -1:
            bounds = self.conf["DoErrorEstimation"]

        # generate states if needed
        if bounds > len(self.mont_carlo_states) - 1:
            self.tomography_states_generator(bounds - len(self.mont_carlo_states) + 1)
        states = np.array(self.mont_carlo_states, dtype="O")[: bounds + 1, :]

        if bounds > 0:
            return getBellSettings_helper_bounds(states[:, 0], rho, partsize_init, partsize, t, bounds)
        else:
            return getBellSettings_helper(rho, partsize_init, partsize, t)

    """
    printLastOutput(bounds)
    Desc: Prints the properties of the last tomography to the console. Properties are defined in tomography conf settings. The calculated properties are determined by self.err_functions.

    Parameters
    ----------
    bounds : int (optional)
        The number of monte carlo runs you want to perform to get a better estimate of each property. Default will use whatever is set in the conf settings.
    """

    def printLastOutput(self, bounds=-1):
        p = np.array(self.last_rho.copy(), dtype="O")
        print("State: ")
        mx = 0
        for i in range(p.shape[0]):
            for j in range(p.shape[1]):
                p[i, j] = floatToString(p[i, j]).replace(" ", "") + "  "
                if len(p[i, j]) > mx:
                    mx = len(p[i, j])

        for i in range(p.shape[0]):
            for j in range(p.shape[1]):
                print(p[i, j] + " " * (mx - len(p[i, j])), end="")
            print("")

        # print(p)
        properties = self.getProperties(bounds)
        for prop in properties:
            if (len(prop) >= 3) and prop[2] != "NA":
                print(prop[0] + " : " + floatToString(prop[1]) + " +/- " + floatToString(prop[2]))
            else:
                print(prop[0] + " : " + floatToString(prop[1]))

    """
    exportToEval(filePath)
    Desc: Exports the input data to the specified file path. You can rerun tomography on this file using importEval() command.

    Parameters
    ----------
    filePath : str (optional)
        The file name you want the data saved to.

    See Also
     ------ 
        importEval
    """

    def exportToEval(self, filePath="pythonEval.txt"):
        TORREPLACE = ""
        # Conf settings
        for k in OLD_FORMAT_CONFIG_KEYS:
            print(k)
            print(isinstance(self.conf[k],np.ndarray))
            print(type(self.conf[k]))
            print(self.conf[k])
            if k == "method":
                TORREPLACE += f'conf["{str(k)}"] = "{str(self.conf[k])}"\n'
            elif isinstance(self.conf[k], np.ndarray):
                A = self.conf[k]
                TORREPLACE += "conf['" + str(k) + "'] = np.array(["
                for i in range(A.shape[0]):
                    if len(A.shape) == 2:
                        TORREPLACE += "["
                        for j in range(A.shape[1]):
                            TORREPLACE += str(A[i, j]) + ","
                        TORREPLACE = TORREPLACE[:-1] + "],"
                    else:
                        TORREPLACE += str(A[i]) + ","
                if len(A.shape) == 2:
                    TORREPLACE = TORREPLACE[:-1] + "])\n"
                else:
                    TORREPLACE = TORREPLACE[:-1] + "])\n"
            else:
                TORREPLACE += f"conf['{str(k)}'] = {str(self.conf[k])}\n"

        # Tomoinput
        A = self.last_input
        TORREPLACE += "tomo_input = np.array(["
        for i in range(A.shape[0]):
            TORREPLACE += "["
            for j in range(A.shape[1]):
                TORREPLACE += str(A[i, j]) + ","
            TORREPLACE = TORREPLACE[:-1] + "],"

        TORREPLACE = TORREPLACE[:-1] + "])\n"
        # intensity
        A = self.intensities
        TORREPLACE += "intensity = np.array(["
        for i in range(A.shape[0]):
            TORREPLACE += str(A[i]) + ","
        TORREPLACE = TORREPLACE[:-1] + "])"

        # print contents to file
        with open(filePath, "w") as f:
            f.write(TORREPLACE)

    """
    export_to_conf_toml(filePath)
    Desc: Exports the configuration that was used by the tomography.
        This exports it in the new format (toml).

    Parameters
    ----------
    filePath : str (optional)
        The file name you want the data saved to.

    See Also
     ------ 
        exportToConf_old_format;import_conf;
    """

    def export_to_conf_toml(self, filePath="pythonConf.toml"):
        new_format_dict = {}

        # Check to see if values exist, otherwise the custom ConfDict
        # will complain about None types being passed in
        for key, value in NEW_FORMAT_CONFIG_KEY_MAPPING.items():
            if isinstance(value, list):
                for val in value:
                    new_format_dict[val] = self.conf[val]
            elif value is not None:
                if self.conf[value] is not None:
                    new_format_dict[key] = self.conf[value]

        # print contents to file
        with open(filePath, "wb") as f:
            tomli_w.dump(new_format_dict, f)

    ## In progress data json construction
    # def export_to_data_json(self, filePath="pythonData.json"):
    #     data_dict = {}
    #     for key, value in NEW_FORMAT_DATA_KEY_MAPPING:
    #         data_dict[key] = self.conf[value]

    #     [coincidences, measurements_densities, measurements_pures, accidentals,overall_norms] = self.filter_data(self.tomo_input)
    #     measurement_states = {}
    #     for i, measurement in enumerate(measurements_pures):
    #         measurement_states[str(i)] = measurement.tolist()

    #     data_dict = {}

    #     with open(filePath, 'w') as f:

    """
    exportToConf(filePath)
    Desc: Exports the conf data to the specified file path. You can rerun tomography on this file using importConf() and importData() command.
        Note: This is being deprecated. You should use the new file format (toml) for configurations.
    Parameters
    ----------
    filePath : str (optional)
        The file name you want the data saved to.

    See Also
     ------ 
        exportToData;importConf;importData
    """

    def exportToConf(self, filePath="pythonConf.txt"):
        TORREPLACE = ""
        # Conf settings
        for k in OLD_FORMAT_CONFIG_KEYS:
            if k == "method":
                TORREPLACE += f'conf["{str(k)}"] = "{str(self.conf[k])}"\n'
            elif isinstance(self.conf[k], np.ndarray):
                A = self.conf[k]
                TORREPLACE += "conf['" + str(k) + "'] = np.array(["
                for i in range(A.shape[0]):
                    if len(A.shape) == 2:
                        TORREPLACE += "["
                        for j in range(A.shape[1]):
                            TORREPLACE += str(A[i, j]) + ","
                        TORREPLACE = TORREPLACE[:-1] + "],"
                    else:
                        TORREPLACE += str(A[i]) + ","
                if len(A.shape) == 2:
                    TORREPLACE = TORREPLACE[:-1] + "])\n"
                else:
                    TORREPLACE = TORREPLACE[:-1] + "])\n"
            else:
                TORREPLACE += f"conf['{str(k)}'] = {str(self.conf[k])}\n"

            # print contents to file
            with open(filePath, "w") as f:
                f.write(TORREPLACE)

    """
     exportToData(filePath)
     Desc: Exports the tomo_input matrix  data to the specified file path. You can rerun tomography on this file using importConf() and importData() command.
        Note: This is being deprecated. You should use the new file format (json) for data.
     Parameters
     ----------
     filePath : str (optional)
         The file name you want the data saved to.

     See Also
      ------ 
         exportToConf;importConf;importData
     """

    def exportToData(self, filePath="pythonData.txt"):
        TORREPLACE = ""

        # Tomoinput
        A = self.last_input
        TORREPLACE += "tomo_input = np.array(["
        for i in range(A.shape[0]):
            TORREPLACE += "["
            for j in range(A.shape[1]):
                TORREPLACE += str(A[i, j]) + ","
            TORREPLACE = TORREPLACE[:-1] + "],"

        TORREPLACE = TORREPLACE[:-1] + "])\n"
        # intensity
        A = self.intensities
        TORREPLACE += "intensity = np.array(["
        for i in range(A.shape[0]):
            TORREPLACE += str(A[i]) + ","
        TORREPLACE = TORREPLACE[:-1] + "])"

        # print contents to file
        with open(filePath, "w") as f:
            f.write(TORREPLACE)

    """
    exportToConf_web(filePath)
    Desc: Exports the conf data to the specified file path. 
    You can run this tomography on our tomography interface at http://tomography.web.engr.illinois.edu/TomographyDemo.php
    THIS IS A WORK IN PROGRESS.

    Parameters
    ----------
    filePath : str (optional)
        The file name you want the data saved to.

    See Also
     ------ 
        exportToData_web
    """

    def exportToConf_web(self, filePath="Config_web.txt"):
        # Conf settings
        TORREPLACE = (
            "conf.NQubits=" + str(self.conf["NQubits"]) + ";\n"
            "conf.NDetectors=" + str(self.conf["NDetectors"]) + ";\n"
            "conf.DoErrorEstimation=" + str(self.conf["DoErrorEstimation"]) + ";\n"
            "conf.Crosstalk=["
        )
        print(self.conf["crosstalk"])
        A = self.conf["crosstalk"]
        for i in range(len(A)):
            TORREPLACE += "["
            for j in range(len(A[i])):
                TORREPLACE += str(A[i][j]).replace(" ", "") + ","
            TORREPLACE = TORREPLACE[:-1] + "],"
        TORREPLACE = TORREPLACE[:-1] + "];\n"

        TORREPLACE += "conf.UseDerivative=0;\n"
        if self.conf["Bellstate"]:
            TORREPLACE += "conf.Bellstate=1;\n"
        else:
            TORREPLACE += "conf.Bellstate=0;\n"
        A = self.conf["Window"]
        if not (isinstance(A, np.ndarray) or isinstance(A, list)):
            A = [A]
        TORREPLACE += "conf.Window=["
        for i in range(len(A)):
            TORREPLACE += str(A[i]).replace(" ", "") + ","
        TORREPLACE = TORREPLACE[:-1] + "];\n"
        A = self.conf["Efficiency"]
        if not (isinstance(A, np.ndarray) or isinstance(A, list)):
            A = [A]
        TORREPLACE += "conf.Efficiency=["
        for i in range(len(A)):
            TORREPLACE += str(A[i]).replace(" ", "") + ","
        TORREPLACE = TORREPLACE[:-1] + "];\n"

        # print(TORREPLACE)
        # print contents to file
        with open(filePath, "w") as f:
            f.write(TORREPLACE)

    """
    exportToData_web(filePath)
    Desc: Exports the tomo_input matrix data to the specified file path. 
    You can run this tomography on our tomography interface at http://tomography.web.engr.illinois.edu/TomographyDemo.php
    THIS IS A WORK IN PROGRESS.
    
    Parameters
    ----------
    filePath : str (optional)
        The file name you want the data saved to.

    See Also
     ------ 
        exportToConf_web
    """

    def exportToData_web(self, filePath="pythonData.txt"):
        TORREPLACE = ""

        # Tomoinput
        A = self.last_input
        TORREPLACE += "tomo_input=["
        for i in range(A.shape[0]):
            TORREPLACE += "["
            for j in range(A.shape[1]):
                element = A[i, j]
                if isinstance(element, complex) and element.imag == 0:
                    element = element.real
                    if modf(element)[0] == 0:
                        element = int(element)

                TORREPLACE += str(element) + ","
            TORREPLACE = TORREPLACE[:-1] + "],"

        TORREPLACE = TORREPLACE[:-1] + "];\n"
        # intensity
        A = self.intensities
        TORREPLACE += "intensity=["
        for i in range(A.shape[0]):
            TORREPLACE += str(A[i]) + ","
        TORREPLACE = TORREPLACE[:-1] + "];"

        # print contents to file
        with open(filePath, "w") as f:
            f.write(TORREPLACE)

    # def exportToData_web(self, filePath="pythonData.txt"):
    #     TORREPLACE = ""
    #
    #     # Tomoinput
    #     A = self.last_input.astype("O")
    #     n_qubit = self.conf['NQubits']
    #     if self.conf['NDetectors'] == 2:
    #         # tomo_input[ :, 1 : 2 * n_qubit + 1 ]: singles
    #         A[:, 1: 2 * n_qubit + 1] = self.last_input[:, 1: 2 * n_qubit + 1].real.astype(int)
    #         # tomo_input[ :, 2 * n_qubit + 1 : 2 ** n_qubit + 2 * n_qubit + 1 ]: coincidences
    #         A[:, 2 * n_qubit + 1: 2 ** n_qubit + 2 * n_qubit + 1] = self.last_input[:,
    #                                                                 2 * n_qubit + 1: 2 ** n_qubit + 2 * n_qubit + 1].real.astype(
    #             int)
    #     else:
    #         # tomo_input[:, 1: n_qubit + 1]: singles
    #         A[:, 1: n_qubit + 1] = self.last_input[:, 1: n_qubit + 1].real.astype(int)
    #         # tomo_input[:, n_qubit + 1]: coincidences
    #         A[:, n_qubit + 1] = self.last_input[:, n_qubit + 1].real.astype(int)
    #
    #     TORREPLACE += "tomo_input=["
    #     for i in range(A.shape[0]):
    #         TORREPLACE += "["
    #         for j in range(A.shape[1]):
    #             if isinstance(A[i, j],int):
    #                 TORREPLACE += str(A[i, j]) + ","
    #             else:
    #                 TORREPLACE += floatToString(A[i, j]) + ","
    #         TORREPLACE = TORREPLACE[:-1] + "],"
    #
    #     TORREPLACE = TORREPLACE[:-1] + "];\n"
    #
    #     # intensity
    #     A = self.intensities
    #     TORREPLACE += "intensity=["
    #     for i in range(A.shape[0]):
    #         TORREPLACE += str(A[i]) + ","
    #     TORREPLACE = TORREPLACE[:-1] + "];"
    #
    #     # print contents to file
    #     with open(filePath, 'w') as f:
    #         f.write(TORREPLACE)
