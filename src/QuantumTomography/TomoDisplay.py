from __future__ import print_function
import scipy as sp
from numpy.core.defchararray import add
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
from mpl_toolkits.mplot3d import Axes3D
from .TomoDisplayHelpers import *
from .TomoFunctions import removeGlobalPhase
from matplotlib.colors import LinearSegmentedColormap
import warnings

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


"""
makeRhoImages(p, plt_given, customColor)
Desc: Creates matlab plots of the density matrix.

Parameters
----------
p : ndarray with shape = (n, 2^numQubits, 2^numQubits)
    The density matrix you want to create plots of.
plt_given : matplotlib.pyplot
    Input pyplot for which the figures will be saved on to.
customColor : boolean
    Specify if you want our custom colorMap. Default is true

See Also
 ------ 
saveRhoImages
"""
def makeRhoImages(p, plt_given, customColor = True):
    # Set up
    numQubits = int(np.log2(p.shape[0]))
    xpos = np.zeros_like(p.flatten(), dtype = float)
    ypos = np.zeros_like(p.flatten(), dtype = float)
    for i in range(0, 2**numQubits):
        xpos[i*2**numQubits:(1+i)*2**numQubits] = .5+i
    for i in range(0, 2**numQubits):
        ypos[i::2**numQubits] = .5+i
    zpos = np.zeros_like(p.flatten(), dtype = float)
    # width of cols
    dx = .9*np.ones_like(xpos)
    dy = .9*np.ones_like(ypos)
    # custom color map
    n_bin = 100
    if(customColor):
        from matplotlib.colors import LinearSegmentedColormap
        cmap_name = 'my_list'
        colors = [(1 / 255.0, 221 / 255.0, 137 / 255.0),
                  (32 / 255.0, 151 / 255.0, 138 / 255.0),
                  (53 / 255.0, 106 / 255.0, 138 / 255.0),
                  (86 / 255.0, 33 / 255.0, 139 / 255.0),
                  (131 / 255.0, 75 / 255.0, 114 / 255.0),
                  (173 / 255.0, 114 / 255.0, 90 / 255.0),
                  (253 / 255.0, 187 / 255.0, 45 / 255.0)]
        colorMap = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)
    else:
        colorMap = plt.cm.jet
    norm = mpl.colors.Normalize(vmin = -1, vmax = 1)

    tickBase = ["H", "V"]
    tick = [""]
    for x in range(numQubits):
        newTick = np.zeros(len(tick)*2, dtype = "O")
        for i in range(len(tick)):
            for j in range(len(tickBase)):
                newTick[len(tick)*i +j] = tick[i] + tickBase[j]
        tick = newTick
    xTicks = ["|"+x+">" for x in tick]
    yTicks = ["|"+x+">" for x in tick]


    # Real Graph
    fig = plt_given.figure()
    ax1 = fig.add_subplot(111, projection = '3d')
    dz = p.flatten().astype(float)
    img = ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color = colorMap((dz + 1) / 2), edgecolor = "black", alpha = .8)



    ax1.axes.set_xticklabels(xTicks)
    ax1.axes.set_yticklabels(yTicks)
    ax1.axes.set_xticks(range(1, 2**numQubits+1))
    ax1.axes.set_yticks(range(1, 2**numQubits+1))
    ax1.axes.set_zticks(np.arange(-1, 1.1, .2))
    ax1.axes.set_zlim3d(-1, 1)
    plt_given.title("Rho Real")
    fig.subplots_adjust(bottom = 0.2)
    ax1 = fig.add_axes([0.2, 0.10, 0.7, 0.065])
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap = colorMap,
                                    norm = norm,
                                    orientation = 'horizontal')

    # Imaginary graph
    fig = plt_given.figure()
    ax1 = fig.add_subplot(111, projection = '3d')
    dz = p.flatten().imag.astype(float)
    ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color = colorMap((dz + 1) / 2), edgecolor = "black", alpha = .8)

    ax1.axes.set_xticklabels(xTicks)
    ax1.axes.set_yticklabels(yTicks)
    ax1.axes.set_xticks(range(1, 2 ** numQubits + 1))
    ax1.axes.set_yticks(range(1, 2 ** numQubits + 1))
    ax1.axes.set_zticks(np.arange(-1, 1.1, .2))
    ax1.axes.set_zlim3d(-1, 1)
    plt_given.title("Rho Imaginary")

    fig.subplots_adjust(bottom = 0.2)
    ax1 = fig.add_axes([0.2, 0.10, 0.7, 0.065])
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap = colorMap,
                                    norm = norm,
                                    orientation = 'horizontal')

"""
saveRhoImages(p, pathToDirectory, customColor)
Desc: Creates and saves matlab plots of the density matrix.

Parameters
----------
p : ndarray with shape = (n, 2^numQubits, 2^numQubits)
    The density matrix you want to create plots of.
pathToDirectory : string
    Path to where you want your images to be saved.

See Also
 ------ 
makeRhoImages
"""
def saveRhoImages(p, pathToDirectory):
    # Set up
    numQubits = int(np.log2(p.shape[0]))
    xpos = np.zeros_like(p.flatten(), dtype = float)
    ypos = np.zeros_like(p.flatten(), dtype = float)
    for i in range(0, 2 ** numQubits):
        xpos[i * 2 ** numQubits:(1 + i) * 2 ** numQubits] = .5 + i
    for i in range(0, 2 ** numQubits):
        ypos[i::2 ** numQubits] = .5 + i
    zpos = np.zeros_like(p.flatten(), dtype = float)
    # width of cols
    dx = .9 * np.ones_like(xpos)
    dy = .9 * np.ones_like(ypos)
    # custom color map
    n_bin = 100
    cmap_name = 'my_list'
    colors = [(1 / 255.0, 221 / 255.0, 137 / 255.0),
              (32 / 255.0, 151 / 255.0, 138 / 255.0),
              (53 / 255.0, 106 / 255.0, 138 / 255.0),
              (86 / 255.0, 33 / 255.0, 139 / 255.0),
              (131 / 255.0, 75 / 255.0, 114 / 255.0),
              (173 / 255.0, 114 / 255.0, 90 / 255.0),
              (253 / 255.0, 187 / 255.0, 45 / 255.0)]
    colorMap = LinearSegmentedColormap.from_list(cmap_name, colors, N = n_bin)

    norm = mpl.colors.Normalize(vmin = -1, vmax = 1)

    tickBase = ["H", "V"]
    tick = [""]
    for x in range(numQubits):
        newTick = np.zeros(len(tick) * 2, dtype = "O")
        for i in range(len(tick)):
            for j in range(len(tickBase)):
                newTick[len(tick) * i + j] = tick[i] + tickBase[j]
        tick = newTick
    xTicks = ["|" + x + ">" for x in tick]
    yTicks = ["|" + x + ">" for x in tick]

    # Real Graph
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection = '3d')
    dz = p.flatten().astype(float)
    img = ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color = colorMap((dz + 1) / 2), edgecolor = "black", alpha = .8)

    ax1.axes.set_xticklabels(xTicks)
    ax1.axes.set_yticklabels(yTicks)
    ax1.axes.set_xticks(range(1, 2 ** numQubits + 1))
    ax1.axes.set_yticks(range(1, 2 ** numQubits + 1))
    ax1.axes.set_zticks(np.arange(-1, 1.1, .2))
    ax1.axes.set_zlim3d(-1, 1)
    plt.title("Rho Real")
    fig.subplots_adjust(bottom = 0.2)
    ax1 = fig.add_axes([0.2, 0.10, 0.7, 0.065])
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap = colorMap,
                                    norm = norm,
                                    orientation = 'horizontal')
    plt.savefig(pathToDirectory + "/rhobarReal.png", bbox_inches = 'tight', pad_inches = 0)

    # Imaginary graph
    fig = plt.figure()
    ax1 = fig.add_subplot(111, projection = '3d')
    dz = p.flatten().imag.astype(float)
    ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color = colorMap((dz + 1) / 2), edgecolor = "black", alpha = .8)

    ax1.axes.set_xticklabels(xTicks)
    ax1.axes.set_yticklabels(yTicks)
    ax1.axes.set_xticks(range(1, 2 ** numQubits + 1))
    ax1.axes.set_yticks(range(1, 2 ** numQubits + 1))
    ax1.axes.set_zticks(np.arange(-1, 1.1, .2))
    ax1.axes.set_zlim3d(-1, 1)
    plt.title("Rho Imaginary")
    fig.subplots_adjust(bottom = 0.2)
    ax1 = fig.add_axes([0.2, 0.10, 0.7, 0.065])
    cb1 = mpl.colorbar.ColorbarBase(ax1, cmap = colorMap,
                                    norm = norm,
                                    orientation = 'horizontal')

    plt.savefig(pathToDirectory + "/rhobarImag.png", bbox_inches = 'tight', pad_inches = 0)


"""
printLastOutput(tomo,bounds)
Desc: Prints the properties of the last tomography to the console. Properties are defined in tomography conf settings. The calculated properties are determined by self.err_functions.

Parameters
----------
tomo : Tomography
    The tomography object you want to see the output of.
bounds : int (optional)
    The number of monte carlo runs you want to perform to get a better estimate of each property. Default will use whatever is set in the conf settings.
"""
def printLastOutput(tomo, bounds = -1):
    p = np.array(tomo.last_rho.copy(), dtype = "O")
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
    properties = tomo.getProperties(bounds)
    for prop in properties:
        if(len(prop) >=3) and prop[2] != 'NA':
            print(prop[0] + " : " + floatToString(prop[1]) + " +/- " + floatToString(prop[2]))
        else:
            print(prop[0] + " : " + floatToString(prop[1]))

"""
matrixToHTML(M)
Desc: Creates an HTML table based on the given matrix.

Parameters
----------
M : 2d numpy array with shape = (2^numQubits, 2^numQubits)
    Matrix you would like to display on your html page.
printEigenVals : boolean
    Specify if you want eigen values to be calculated and displayed at the bottom of the table.
Returns
-------
res : string
    HTML code of the created table.

See Also
 ------ 
propertiesToHTML
"""
def matrixToHTML(M, printEigenVals = False):
    s = np.shape(M)
    res = '<table class="KwiatDataMatrix" style = \"border: 1px solid black;border-collapse: collapse;font-size: 15px; table-layout:fixed;width:100%;margin-top: 25px;\">'
    for i in range(s[0]):
        res = res+' <tr>'
        for j in range(s[1]):
            # res = res + '<td style = "border: 1px solid black;">' + str(np.real(M[i, j])) + "<div style = \"color:rebeccapurple;font-weight: bold;display:inline;\">+</div><BR>"+ str(np.imag(M[i, j]))
            # res = res + '<div style = \"color:rebeccapurple;font-weight: bold;display:inline;\">j</div></td>'
            res = res + '<td style = "border: 1px solid black;">' + floatToString(M[i, j], True)+ '</td>'
        res = res +'</tr>'
    res = res+'</table>'
    if(printEigenVals):
        d, v = np.linalg.eig(M)
        sum = 0
        eigenVals = "<h5>Eigen Values : "
        for x in range(0, len(d)):
            eigenVals = eigenVals+str(round(d[x].real, 5))
            if(abs(d[x].imag)>.00001):
                eigenVals = eigenVals+"< div style = \"color:rebeccapurple;font-weight: bold;display:inline;\">+</div>"
                eigenVals = eigenVals+str(round(d[x].imag, 5))
                eigenVals = eigenVals+"<div style = \"color:rebeccapurple;font-weight: bold;display:inline;\">j</div>"
            eigenVals = eigenVals+" , "
            sum+= d[x]
        eigenVals = str(eigenVals)[0:len(str(eigenVals))-2]
        eigenVals = eigenVals +"</h5>"
        res = res+eigenVals
    return res


"""
propertiesToHTML(vals)
Desc: Creates an HTML table based on the given property values.

Parameters
----------
vals : ndarray with shape = (length of self.err_functions, 2)
    The first col is the name of the property.
    The second col is the value of the property.
    The third col is the error bound on the property.

Returns
-------
res : string
    HTML code of the created table.

See Also
 ------ 
matrixToHTML
"""
def propertiesToHTML(vals):
    f = '<h3 >Properties of Rho</h3><table style = \"width:60%;margin-top:10px;font-size: 15px;padding-bottom:5px;float:none;\"><tr><td style = "font-size: 20;font-weight: 1000;color: rebeccapurple;">Property</td>'
    hasSTD = len(vals) > 2
    if hasSTD and any(vals[:,2] != "NA"):
        f += '<td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">Average Value</td>'
        f += '<td style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">STD Error</td></tr>'
        f += '<td style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;"></td></tr>'
    else:
        f += '<td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">Value</td>'
        f += '<td style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;"></td></tr>'
    for v in vals:
        if(v[1] != "NA"):
            f += '<tr>'
            f += '<td><div onmouseover = "Tip( ' + v[0].replace(" ", "") + 'Tip)" onmouseout = "hideTip()">'+v[0] + '</td>'
            f += '<td name = "' + v[0].replace(" ", "") + '_value">' + floatToString(v[1], True) + '</td>'
            if hasSTD and v[2] != "NA":
                f += '<td> +/- ' + floatToString(v[2], True) + '</td>'
            else:
                f += "<td></td>"
            f += '</tr>'
    f += "</table>"
    return f


"""
stateToString(vals)
Desc: Creates a string of the 1d pure state.

Parameters
----------    
pure_state : 1darray with length = 2
    The state in ket form.
"""
def stateToString(state):
    state = removeGlobalPhase(state)
    return floatToString(state[0])+"|H> + "+floatToString(state[1])+"|V>"
