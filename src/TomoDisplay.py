import matplotlib.pyplot as plt
import numpy as np

def showRhoImages(p):
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    xpos = [.5, .5, .5, .5,
            1.5, 1.5, 1.5, 1.5,
            2.5, 2.5, 2.5, 2.5,
            3.5, 3.5, 3.5, 3.5]
    ypos = [.5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5, ]
    zpos = np.zeros(16)
    dx = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]
    dy = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]

    dz = [.5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5, ]
    dz = p.flatten().astype(float)


    colors = ['r','g','b','purple']

    #ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
    for i in range(4):
        ax1.bar3d(xpos[i*4:i*4+4], ypos[i*4:i*4+4], zpos[i*4:i*4+4], dx[i*4:i*4+4], dy[i*4:i*4+4], dz[i*4:i*4+4],edgecolor="black",alpha=.65, color=colors[i])

    ax1.axes.set_xticklabels(["|HH>", "|HV>", "|VH>", "|VV>"])
    ax1.axes.set_yticklabels(["<HH|", "<HV|", "<VH|", "<VV|"])
    ax1.axes.set_xticks([1, 2, 3, 4])
    ax1.axes.set_yticks([1, 2, 3, 4])
    ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
    ax1.axes.set_zlim3d(-1, 1)
    plt.title("Rho")



    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection='3d')
    xpos = [.5, .5, .5, .5,
            1.5, 1.5, 1.5, 1.5,
            2.5, 2.5, 2.5, 2.5,
            3.5, 3.5, 3.5, 3.5]
    ypos = [.5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5,
            .5, 1.5, 2.5, 3.5, ]
    zpos = np.zeros(16)
    dx = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]
    dy = [.9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9
        , .9, .9, .9, .9]

    dz = [.5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5,
          .5, 1.5, 2.5, 3.5, ]
    dz = p.flatten().imag.astype(float)
    colors = ['r', 'g', 'b', 'purple']

    # ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
    for i in range(4):
        ax1.bar3d(xpos[i * 4:i * 4 + 4], ypos[i * 4:i * 4 + 4], zpos[i * 4:i * 4 + 4], dx[i * 4:i * 4 + 4],
                  dy[i * 4:i * 4 + 4], dz[i * 4:i * 4 + 4], edgecolor="black", alpha=.65, color=colors[i])


    ax1.axes.set_xticklabels(["|HH>", "|HV>", "|VH>", "|VV>"])
    ax1.axes.set_yticklabels(["<HH|", "<HV|", "<VH|", "<VV|"])
    ax1.axes.set_xticks([1, 2, 3, 4])
    ax1.axes.set_yticks([1, 2, 3, 4])
    ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
    ax1.axes.set_zlim3d(-1, 1)
    plt.title("Rho Imaginary")
    plt.show()

def displayOutput(p,inten,fval,tomo):
    err_n = tomo.conf['DoErrorEstimation']
    vals = tomo.getproperties(p)
    if (err_n > 1):
        rhon = tomo.tomography_error_states_generator(err_n)
        [mean, errs, mean_fid, err_fid] = tomo.tomography_error(p, rhon)
        errs = np.around(errs, 5)
        vals = np.around(vals, 5)
        intensity = round(inten, 1)
        fval = round(fval, 5)
        outputVals, outputNames = outputvalues(tomo.err_functions, vals, errs, intensity, fval, err_n)
    else:
        outputNames = tomo.err_functions
        outputVals = {}
        for x in range(0, len(outputNames)):
            outputVals[outputNames[x]] = [vals[x]]

    for x in range(0, len(tomo.err_functions)):
        if(len(outputVals[outputNames[x]]) == 1):
            print(outputNames[x], " = ", outputVals[outputNames[x]][0])
        else:
            print(outputNames[x], " = ", outputVals[outputNames[x]][0], "+/-", outputVals[outputNames[x]][1])
    print("State = ")
    print(p)
    return outputVals, outputNames
def outputvalues(func, values, errors, inten, fvalp, err_time):
    output = {"Intensity": (inten,), "fval": (fvalp,), "Error Estimation Times": (err_time,)}
    name = ["Intensity", "fval", "Error Estimation Times"]
    for i in range(len(func)):
        if func[i] == 'concurrence':
            output["Concurrence"] = (values[i], errors[i])
            name = np.append(name, "Concurrence")
        elif func[i] == 'tangle':
            output["Tangle"] = (values[i], errors[i])
            name = np.append(name, "Tangle")
        elif func[i] == 'entanglement':
            output["Entanglement"] = (values[i], errors[i])
            name = np.append(name, "Entanglement")
        elif func[i] == 'entropy':
            output["Entropy"] = (values[i], errors[i])
            name = np.append(name, "Entropy")
        elif func[i] == 'linear_entropy':
            output["Linear Entropy"] = (values[i], errors[i])
            name = np.append(name, "Linear Entropy")
        elif func[i] == 'negativity':
            output["Negativity"] = (values[i], errors[i])
            name = np.append(name, "Negativity")
    return output, name