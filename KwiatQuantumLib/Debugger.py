
# ##################
# ## TESTING DATA ##
# ##################
# #Uncomment from one green line to another to use one section of sample data for testing.
# # Only use one at a time, the program will run everything uncommented but will only use the last section
#
# #Input variables names: 'tomo_input'(2-d array),'intensity'(1-d array),'conf'(dictionary), 'err_n'(int),
# # 'err_functions'(string array).
#
# #More details about each input can be found in the README pdf tutorial on our website
#
# """SAMPLE SET 1"""
# # #Basic settings and Quantum State is |HH>
# # tomo_input = np.array([[1,0,0,500,1,0,1,0],
# #                        [1,0,0,0,1,0,0,1],
# #                        [1,0,0,250,1,0,0.7071,0.7071],
# #                        [1,0,0,250,1,0,0.7071,0.7071j],
# #                        [1,0,0,0,0,1,1,0],[1,0,0,0,0,1,0,1],
# #                        [1,0,0,0,0,1,0.7071,0.7071],
# #                        [1,0,0,0,0,1,0.7071,0.7071j],
# #                        [1,0,0,250,0.7071,0.7071,1,0],
# #                        [1,0,0,0,0.7071,0.7071,0,1],
# #                        [1,0,0,125,0.7071,0.7071,0.7071,0.7071],
# #                        [1,0,0,125,0.7071,0.7071,0.7071,0.7071j],
# #                        [1,0,0,250,0.7071,0.7071j,1,0],
# #                        [1,0,0,0,0.7071,0.7071j,0,1],
# #                        [1,0,0,125,0.7071,0.7071j,0.7071,0.7071],
# #                        [1,0,0,125,0.7071,0.7071j,0.7071,0.7071j]])
# # intensity = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
# # conf = {'NQubits': 2,
# #         'NDetectors': 1,
# #         'Crosstalk': np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]),
# #         'UseDerivative': 0,
# #         'DoErrorEstimation': 0,
# #         'DoDriftCorrection': 'no',
# #         'DoAccidentalCorrection' : 1,
# #         'Window': 0,
# #         'Efficiency': [1, 1, 1, 1]}
#
# """SAMPLE SET 2"""
# # #Quantum State is Bellstate: .7071|HH> + .7071|VV>
# # #Settings: 32 measurements, 2 det/qubit
# # tomo_input = np.array([[1,0,0,0,0,250,0,0,250,1,0,1,0],
# #                         [1,0,0,0,0,0,250,250,0,1,0,0,1],
# #                         [1,0,0,0,0,125,125,125,125,1,0,0.7071,0.7071],
# #                         [1,0,0,0,0,125,125,125,125,1,0,0.7071,-0.7071],
# #                         [1,0,0,0,0,125,125,125,125,1,0,0.7071,0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,1,0,0.7071,-0.7071j],
# #                         [1,0,0,0,0,0,250,250,0,0,1,1,0],
# #                         [1,0,0,0,0,250,0,0,250,0,1,0,1],
# #                         [1,0,0,0,0,125,125,125,125,0,1,0.7071,0.7071],
# #                         [1,0,0,0,0,125,125,125,125,0,1,0.7071,-0.7071],
# #                         [1,0,0,0,0,125,125,125,125,0,1,0.7071,0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,0,1,0.7071,-0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071,1,0],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071,0,1],
# #                         [1,0,0,0,0,250,0,0,250,0.7071,0.7071,0.7071,0.7071],
# #                         [1,0,0,0,0,0,250,250,0,0.7071,0.7071,0.7071,-0.7071],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071,0.7071,0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071,0.7071,-0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071,1,0],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071,0,1],
# #                         [1,0,0,0,0,0,250,250,0,0.7071,-0.7071,0.7071,0.7071],
# #                         [1,0,0,0,0,250,0,0,250,0.7071,-0.7071,0.7071,-0.7071],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071,0.7071,0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071,0.7071,-0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071j,1,0],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071j,0,1],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071j,0.7071,0.7071],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,0.7071j,0.7071,-0.7071],
# #                         [1,0,0,0,0,0,250,250,0,0.7071,0.7071j,0.7071,0.7071j],
# #                         [1,0,0,0,0,250,0,0,250,0.7071,0.7071j,0.7071,-0.7071j],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071j,1,0],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071j,0,1],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071j,0.7071,0.7071],
# #                         [1,0,0,0,0,125,125,125,125,0.7071,-0.7071j,0.7071,-0.7071],
# #                         [1,0,0,0,0,250,0,0,250,0.7071,-0.7071j,0.7071,0.7071j],
# #                         [1,0,0,0,0,0,250,250,0,0.7071,-0.7071j,0.7071,-0.7071j]])
# # intensity = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
# # conf =  {'NQubits': 2,
# #          'NDetectors': 2,
# #          'Crosstalk': np.array([[1, 0, 0, 0],[0, 1, 0, 0],[0, 0, 1, 0],[0, 0, 0, 1]]),
# #          'UseDerivative': 0,
# #          'DoErrorEstimation': '0',
# #          'DoDriftCorrection': 'no',
# #          'DoAccidentalCorrection' : 1,
# #          'Window': [1, 1, 1, 1],
# #          'Efficiency': [1, 1, 1, 1]}
#
# """DEBUG SET 1"""
# # # Basic settings and random quantum state and create data
# #
# # hh = np.random.randint(-1, 2) * np.random.random() + np.random.randint(-1, 2) * np.random.random() * 1j
# # hv = np.random.randint(-1, 2) * np.random.random() + np.random.randint(-1, 2) * np.random.random() * 1j
# # vh = np.random.randint(-1, 2) * np.random.random() + np.random.randint(-1, 2) * np.random.random() * 1j
# # vv = np.random.randint(-1, 2) * np.random.random() + np.random.randint(-1, 2) * np.random.random() * 1j
# # sum = hh * np.conjugate(hh) + hv * np.conjugate(hv) + vh * np.conjugate(vh) + vv * np.conjugate(vv)
# # hh = hh / np.sqrt(sum)
# # hv = hv / np.sqrt(sum)
# # vh = vh / np.sqrt(sum)
# # vv = vv / np.sqrt(sum)
# #
# # tomo_input = np.array([[1,0,0,0,1,0,1,0],
# #                        [1,0,0,0,1,0,0,1],
# #                        [1,0,0,0,1,0,0.7071,0.7071],
# #                        [1,0,0,0,1,0,0.7071,0.7071j],
# #                        [1,0,0,0,0,1,1,0],[1,0,0,0,0,1,0,1],
# #                        [1,0,0,0,0,1,0.7071,0.7071],
# #                        [1,0,0,0,0,1,0.7071,0.7071j],
# #                        [1,0,0,0,0.7071,0.7071,1,0],
# #                        [1,0,0,0,0.7071,0.7071,0,1],
# #                        [1,0,0,0,0.7071,0.7071,0.7071,0.7071],
# #                        [1,0,0,0,0.7071,0.7071,0.7071,0.7071j],
# #                        [1,0,0,0,0.7071,0.7071j,1,0],
# #                        [1,0,0,0,0.7071,0.7071j,0,1],
# #                        [1,0,0,0,0.7071,0.7071j,0.7071,0.7071],
# #                        [1,0,0,0,0.7071,0.7071j,0.7071,0.7071j]])
# # intensity = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
# # conf = {'NQubits': 2,
# #         'NDetectors': 1,
# #         'Crosstalk': np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]),
# #         'UseDerivative': 0,
# #         'DoErrorEstimation': 0,
# #         'DoDriftCorrection': 'no',
# #          'DoAccidentalCorrection' : 1,
# #         'Window': 0,
# #         'Efficiency': [1, 1, 1, 1]}
# # x = 0
# # while x < len(tomo_input):
# #     tomo_input[x][3] += hh * tomo_input[x][4] * tomo_input[x][6]
# #     tomo_input[x][3] += hv * tomo_input[x][4] * tomo_input[x][7]
# #     tomo_input[x][3] += vh * tomo_input[x][5] * tomo_input[x][6]
# #     tomo_input[x][3] += vv * tomo_input[x][5] * tomo_input[x][7]
# #     tomo_input[x][3] = tomo_input[x][3] * np.conjugate(tomo_input[x][3]) * 500
# #     x += 1
# # # #Calculates the density matrix of the random state for debuggin
# # # temp1 = [[np.conjugate(hh),
# # #           np.conjugate(hv),
# # #           np.conjugate(vh),
# # #           np.conjugate(vv)]]
# # # temp2 = [[hh],
# # #          [hv],
# # #          [vh],
# # #          [vv]
# # #          ]
# # # densityMatrix = tensor_product(temp2, temp1)
# # # #^Keep commented unless debugging
#
# """DEBUG SET 2"""
# # You set quantum state and create data for 2 qubits
# # #Settings: 32 measurements
# # #Quantum state is in the form : hh|HH> + hv|HV> + vh|VH> + vv|VV>
#
# hh = .5
# hv = .5
# vh = .5
# vv = .5
#
# tomo_input = np.array([[1,0,0,0,1,0,1,0],[1,0,0,0,1,0,0,1],[1,0,0,0,1,0,0.7071,0.7071],[1,0,0,0,1,0,0.7071,-0.7071],[1,0,0,0,1,0,0.7071,0.7071j],[1,0,0,0,1,0,0.7071,-0.7071j],[1,0,0,0,0,1,1,0],[1,0,0,0,0,1,0,1],[1,0,0,0,0,1,0.7071,0.7071],[1,0,0,0,0,1,0.7071,-0.7071],[1,0,0,0,0,1,0.7071,0.7071j],[1,0,0,0,0,1,0.7071,-0.7071j],[1,0,0,0,0.7071,0.7071,1,0],[1,0,0,0,0.7071,0.7071,0,1],[1,0,0,0,0.7071,0.7071,0.7071,0.7071],[1,0,0,0,0.7071,0.7071,0.7071,-0.7071],[1,0,0,0,0.7071,0.7071,0.7071,0.7071j],[1,0,0,0,0.7071,0.7071,0.7071,-0.7071j],[1,0,0,0,0.7071,-0.7071,1,0],[1,0,0,0,0.7071,-0.7071,0,1],[1,0,0,0,0.7071,-0.7071,0.7071,0.7071],[1,0,0,0,0.7071,-0.7071,0.7071,-0.7071],[1,0,0,0,0.7071,-0.7071,0.7071,0.7071j],[1,0,0,0,0.7071,-0.7071,0.7071,-0.7071j],[1,0,0,0,0.7071,0.7071j,1,0],[1,0,0,0,0.7071,0.7071j,0,1],[1,0,0,0,0.7071,0.7071j,0.7071,0.7071],[1,0,0,0,0.7071,0.7071j,0.7071,-0.7071],[1,0,0,0,0.7071,0.7071j,0.7071,0.7071j],[1,0,0,0,0.7071,0.7071j,0.7071,-0.7071j],[1,0,0,0,0.7071,-0.7071j,1,0],[1,0,0,0,0.7071,-0.7071j,0,1],[1,0,0,0,0.7071,-0.7071j,0.7071,0.7071],[1,0,0,0,0.7071,-0.7071j,0.7071,-0.7071],[1,0,0,0,0.7071,-0.7071j,0.7071,0.7071j],[1,0,0,0,0.7071,-0.7071j,0.7071,-0.7071j]])
# intensity = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
# conf = {'NQubits': 2,
#         'NDetectors': 1,
#         'Crosstalk': np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]),
#         'UseDerivative': 0,
#         'DoErrorEstimation': 0,
#         'DoDriftCorrection': 'no',
#         'DoAccidentalCorrection' : 1,
#         'Window': 0,
#         'Efficiency': [1, 1, 1, 1]}
# x = 0
# while x < len(tomo_input):
#     tomo_input[x][3] += hh * tomo_input[x][4] * tomo_input[x][6]
#     tomo_input[x][3] += hv * tomo_input[x][4] * tomo_input[x][7]
#     tomo_input[x][3] += vh * tomo_input[x][5] * tomo_input[x][6]
#     tomo_input[x][3] += vv * tomo_input[x][5] * tomo_input[x][7]
#     tomo_input[x][3] = tomo_input[x][3] * np.conjugate(tomo_input[x][3]) * 500
#     x += 1
# """DEBUG SET 3"""
# # # You set quantum state and create data for 1 qubit
# # # #Settings: 32 measurements
# # # #Quantum state is in the form : hh|HH> + hv|HV> + vh|VH> + vv|VV>
# #
# # hh = .5
# # hv = .5
# # vh = .5
# # vv = .5
# #
# # tomo_input = np.array([[1,0,0,0,1,0,1,0],[1,0,0,0,1,0,0,1],[1,0,0,0,1,0,0.7071,0.7071],[1,0,0,0,1,0,0.7071,-0.7071],[1,0,0,0,1,0,0.7071,0.7071j],[1,0,0,0,1,0,0.7071,-0.7071j],[1,0,0,0,0,1,1,0],[1,0,0,0,0,1,0,1],[1,0,0,0,0,1,0.7071,0.7071],[1,0,0,0,0,1,0.7071,-0.7071],[1,0,0,0,0,1,0.7071,0.7071j],[1,0,0,0,0,1,0.7071,-0.7071j],[1,0,0,0,0.7071,0.7071,1,0],[1,0,0,0,0.7071,0.7071,0,1],[1,0,0,0,0.7071,0.7071,0.7071,0.7071],[1,0,0,0,0.7071,0.7071,0.7071,-0.7071],[1,0,0,0,0.7071,0.7071,0.7071,0.7071j],[1,0,0,0,0.7071,0.7071,0.7071,-0.7071j],[1,0,0,0,0.7071,-0.7071,1,0],[1,0,0,0,0.7071,-0.7071,0,1],[1,0,0,0,0.7071,-0.7071,0.7071,0.7071],[1,0,0,0,0.7071,-0.7071,0.7071,-0.7071],[1,0,0,0,0.7071,-0.7071,0.7071,0.7071j],[1,0,0,0,0.7071,-0.7071,0.7071,-0.7071j],[1,0,0,0,0.7071,0.7071j,1,0],[1,0,0,0,0.7071,0.7071j,0,1],[1,0,0,0,0.7071,0.7071j,0.7071,0.7071],[1,0,0,0,0.7071,0.7071j,0.7071,-0.7071],[1,0,0,0,0.7071,0.7071j,0.7071,0.7071j],[1,0,0,0,0.7071,0.7071j,0.7071,-0.7071j],[1,0,0,0,0.7071,-0.7071j,1,0],[1,0,0,0,0.7071,-0.7071j,0,1],[1,0,0,0,0.7071,-0.7071j,0.7071,0.7071],[1,0,0,0,0.7071,-0.7071j,0.7071,-0.7071],[1,0,0,0,0.7071,-0.7071j,0.7071,0.7071j],[1,0,0,0,0.7071,-0.7071j,0.7071,-0.7071j]])
# # intensity = np.array([1,1,1,1,1,1,1,1,1,1,1,1,1,1,1,1])
# # conf = {'NQubits': 2,
# #         'NDetectors': 1,
# #         'Crosstalk': np.array([[1,0,0,0],[0,1,0,0],[0,0,1,0],[0,0,0,1]]),
# #         'UseDerivative': 0,
# #         'DoErrorEstimation': 0,
# #         'DoDriftCorrection': 'no',
# #         'DoAccidentalCorrection' : 1,
# #         'Window': 0,
# #         'Efficiency': [1, 1, 1, 1]}
# # x = 0
# # while x < len(tomo_input):
# #     tomo_input[x][3] += hh * tomo_input[x][4] * tomo_input[x][6]
# #     tomo_input[x][3] += hv * tomo_input[x][4] * tomo_input[x][7]
# #     tomo_input[x][3] += vh * tomo_input[x][5] * tomo_input[x][6]
# #     tomo_input[x][3] += vv * tomo_input[x][5] * tomo_input[x][7]
# #     tomo_input[x][3] = tomo_input[x][3] * np.conjugate(tomo_input[x][3]) * 500
# #     x += 1
# """Get input from file"""
# # conf = {'NQubits': 2, 'Properties': ['concurrence','tangle','entanglement','entropy','linear_entropy','negativity']}
# #
# # 'QuditSizes'
# #
# # # whole function array is: ['concurrence','tangle','entanglement','entropy','linear_entropy','negativity']
# #
# # part_init = 9
# # part = 5
# # bell_times = 3
# # conf['DoAccidentalCorrection'] = 1
# # #execfile('pythoneval.txt')
# # exec(compile(open('Python_Testing\pythoneval.txt', "rb").read(), 'pythoneval.txt', 'exec'))
# """-----------------------"""
#
# err_n = error_n(conf)
#
# err_functions = getfield_default(conf, 'Properties', [])
#
# [rho, intensity, fval] = state_tomography(tomo_input, intensity, conf)
# ###################
# ## MAIN FUNCTION ##
# ###################
#
# """NUMBER OF ERROR ESTIMATES"""
# # default: err_n = 100
# err_n = 2
#
# """CHOOSE ERROR FUNCTIONS"""
# # All the available error functions are: ['concurrence','tangle','entanglement','entropy','linear_entropy','negativity']
# err_functions = ['tangle', 'linear_entropy', 'entropy']
#
# """MAIN FUNCTIONS"""
# [rho, intensity, fval] = state_tomography(tomo_input, intensity, conf)
#
# rhon = tomography_error_states_generator(tomo_input, intensity, conf, err_n)
# [mean, errs, mean_fid, err_fid] = tomography_error(rho,rhon,err_functions, err_n)
#
# """GRAPH MATRIX AND PRINT ERRORS"""
# displayOutput(rho,intensity,fval,errs,mean,mean_fid)
#

#########################
# COINCIDENCE GENERATOR #
#########################


def measure_2qubit(m1, m2):
    meas = [[1, 0], [0, 1], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2)], [np.sqrt(1.0 / 2), -np.sqrt(1.0 / 2)],
            [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2) * 1j], [np.sqrt(1.0 / 2), -np.sqrt(1.0 / 2) * 1j]]
    meas = tensor_product_2(meas[m1], meas[m2])
    meas = np.outer(meas, meas.transpose().conj())
    return meas


def measure_opt(opt):
    if opt == 6:
        meas = [[1, 0], [0, 1], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2)], [np.sqrt(1.0 / 2), -np.sqrt(1.0 / 2)],
                [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2) * 1j], [np.sqrt(1.0 / 2), -np.sqrt(1.0 / 2) * 1j]]
    elif opt == 4:
        meas = [[1, 0], [0, 1], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2)], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2) * 1j]]
    elif opt == 3:
        meas = [[1, 0], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2)], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2) * 1j]]
    else:
        print('Option is invalid.')
        meas = 0
    return meas


def measure_nqubit(meas, n, m):
    rv = 1
    if n > 1:
        for i in range(n):
            rv = tensor_product_2(meas[m[n - 1 - i]], rv)
    elif n == 1:
        m = m[0]
        rv = np.array(meas[m])

    rv = np.outer(rv, rv.transpose().conj())
    return rv


def measure_generator(nbits=2, ndet=1, option='HVDARL'):
    d = 2 ** nbits
    if ndet == 1:
        if option == 'HVDARL':
            opt = 6
        elif option == 'HVDR':
            opt = 4
        elif option == "HDR":
            opt = 3
        else:
            print('Wrong measurement option format.')
            opt = 0

        n = opt ** nbits
        meas = np.zeros([n, d, d]) + 0j
        m = map(int, np.zeros(nbits))
        m_opt = measure_opt(opt)
        for i in range(n):
            for j in range(nbits):
                if ((np.array(m[j:nbits-1]) == np.repeat(opt-1, nbits-1-j)).all()) & (m[nbits-1] == opt):
                    m[j-1] += 1
                    m[j:nbits] = map(int, np.zeros(nbits-j))
                    break
            meas[i] = measure_nqubit(m_opt, nbits, m)
            m[nbits-1] += 1

    elif ndet == 2:
        n = 3 ** nbits
        meas = np.zeros([n, d, d, d]) + 0j
        m = np.zeros(nbits)
        m_opt = measure_opt(6)
        for i in range(n):
            for j in range(nbits):
                if ((np.array(m[j:nbits-1]) == np.repeat(4, nbits-1-j)).all()) & (m[nbits-1] == 6):
                    m[j-1] += 2
                    m[j:nbits] = map(int, np.zeros(nbits-j))
                    break
            mm = np.zeros(nbits)
            for ii in range(d):
                for jj in range(nbits):
                    if ((np.array(mm[jj:nbits-1]) == np.repeat(1, nbits-1-jj)).all()) & (mm[nbits-1] == 2):
                        mm[jj-1] += 1
                        mm[jj:nbits] = map(int, np.zeros(nbits-jj))
                        break
                mmm = mm + m
                mmm = map(int, mmm)
                meas[i, ii] = measure_nqubit(m_opt, nbits, mmm)
                mm[nbits-1] += 1
            m[nbits-1] += 2

    else:
        print('Wrong detector number.')
        meas = 0

    return meas


def coinc_generator(rhoo, meas, n):
    ndet = len(meas.shape) - 2
    if ndet == 1:
        n_m = meas.shape[0]
        coinc = np.zeros(n_m)
        for i in range(n_m):
            coinc[i] = np.real(np.trace(np.dot(meas[i], rhoo)) * n)
            coinc[i] = np.random.poisson(coinc[i])

    elif ndet == 2:
        n_m = [meas.shape[0], meas.shape[1]]
        coinc = np.zeros(n_m)
        for i in range(n_m[0]):
            for j in range(n_m[1]):
                coinc[i, j] = np.real(np.trace(np.dot(meas[i, j], rhoo)) * n)
                coinc[i, j] = np.random.poisson(coinc[i, j])
    else:
        print('Wrong detector number.')
        coinc = 0

    return coinc


def measurement_opt(opt):
    if opt == 6:
        meas = [[1, 0], [0, 1], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2)], [np.sqrt(1.0 / 2), -np.sqrt(1.0 / 2)],
                [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2) * 1j], [np.sqrt(1.0 / 2), -np.sqrt(1.0 / 2) * 1j]]
    elif opt == 4:
        meas = [[1, 0], [0, 1], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2)], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2) * 1j]]
    elif opt == 3:
        meas = [[1, 0], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2)], [np.sqrt(1.0 / 2), np.sqrt(1.0 / 2) * 1j]]
    else:
        print('Option is mot valid.')
        meas = 0
    return meas


def measurement_generator(n_qubit, option='HVDARL'):
    if option == 'HVDARL':
        opt = 6
    elif option == 'HVDR':
        opt = 4
    elif option == "HDR":
        opt = 3
    else:
        print('Wrong measurement option format.')
        opt = 0
    m_opt = measurement_opt(opt)
    d = 2*n_qubit
    x = opt**n_qubit
    meas = np.zeros([x, d])+0j
    m = map(int, np.zeros(n_qubit))
    for i in range(x):
        for j in range(n_qubit):
            if ((np.array(m[j:n_qubit-1]) == np.repeat(opt-1, n_qubit-1-j)).all()) & (m[n_qubit-1] == opt):
                m[j-1] += 1
                m[j:n_qubit] = map(int, np.zeros(n_qubit-j))
                break
        for j in range(n_qubit):
            meas[i, 2*j] = m_opt[m[j]][0]
            meas[i, 2*j+1] = m_opt[m[j]][1]
        m[n_qubit-1] += 1

    return meas


def random_state(n=2):
    d = 2 ** n
    t = np.zeros([d, d]) + 0j
    for i in range(d):
        for j in range(d):
            if i == j:
                t[i, j] = np.random.rand()
            elif i < j:
                t[i, j] = np.random.rand() + np.random.rand() * 1j

    rhoo = (np.dot(t, t.transpose().conj()))
    rhoo /= np.trace(rhoo)
    return rhoo


def randu(d):
    rm = np.random.randn(d, d) + np.random.randn(d, d) * 1j
    ru = sp.linalg.orth(rm)
    return ru


def random_pure_state(n=2):
    d = 2 ** n
    a = np.zeros((d, d))
    a[0, 0] = 1
    ru = randu(d)
    rhoo = np.dot(ru, np.dot(a, ru.transpose().conj()))
    return rhoo


def crosstalk_generator(confp):
    ctalk = getfield_default(confp, 'ctalk', [[1, 0], [0, 1]])
    nbits = confp['NQubits']
    crosstalk = 1
    for n in range(nbits):
        crosstalk = tensor_product_2(crosstalk, ctalk)

    return crosstalk

##################
# DATA GENERATOR #
##################


def data_generator(rhoo, ndet=1, option='HVDARL', n=1000):
    if option == 'HVDARL':
        opt = 6
    elif option == 'HVDR':
        opt = 4
    elif option == "HDR":
        opt = 3
    else:
        print('Wrong measurement option format.')
        opt = 0

    n_qubit = np.int(np.log2(len(rhoo)))

    if ndet == 1:
        d = opt**n_qubit
        x = 2 + 3 * n_qubit
        input = np.zeros([d, x]) + 0j
        input[:, 0] = np.ones(d)
        input[:, n_qubit+1] = coinc_generator(rhoo, measure_generator(n_qubit, ndet, option), n)
        input[:, np.arange(n_qubit+2, len(input[0]))] = measurement_generator(n_qubit, option)

    elif ndet == 2:
        opt = 3
        d = opt**n_qubit
        n_singles = 2*n_qubit
        n_coinc = 2**n_qubit
        x = 1 + n_coinc + n_singles + 2 * n_qubit
        input = np.zeros([d, x]) + 0j
        input[:, 0] = np.ones(d)
        input[:, np.arange(n_singles+1, n_singles+n_coinc+1)] = \
            coinc_generator(rhoo, measure_generator(n_qubit, ndet, option), n)
        input[:, np.arange(n_singles+n_coinc+1, len(input[0]))] = measurement_generator(n_qubit, 'HDR')

    else:
        print('Wrong detector number.')
        input = 0

    return input


def data_collector(coinc, sings, t, confp):
    ndet = confp['NDetectors']
    nbits = confp['NQubits']
    option = confp['Measurements']
    meas = measurement_generator(nbits, option)
    if option == 'HVDARL':
        opt = 6
    elif option == 'HVDR':
        opt = 4
    elif option == "HDR":
        opt = 3
    else:
        print('Wrong measurement option format.')
        opt = 0

    if ndet == 1:
        d = opt**nbits
        x = 2+3*nbits
        input = np.zeros([d, x]) + 0j
        input[:, 0] = t
        input[:, np.arange(1, nbits+1)] = sings
        input[:, nbits+1] = coinc
        input[:, np.arange(nbits+2, len(input[0]))] = meas

    elif ndet == 2:
        opt = 3
        d = opt**nbits
        x = 2**nbits+4*nbits+1
        input = np.zeros([d, x]) + 0j
        input[:, 0] = t
        input[:, np.arange(1, 2*nbits+1)] = sings
        input[:, np.arange(2*nbits+1, 2**nbits+2*nbits+1)] = coinc
        input[:, np.arange(2**nbits+2*nbits+1, len(input[0]))] = meas
    else:
        print('Wrong detector number.')
        input = 0

    return input