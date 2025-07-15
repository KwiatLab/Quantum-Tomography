from __future__ import print_function
import argparse
import os
import sys
import warnings
warnings.filterwarnings('ignore')
import numpy as np
from numpy.core.defchararray import add
import QuantumTomography as qKLib

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""

# This script is used to run quantum tomography from the command line.
#
# Running the command Quantum-Tomography from the command line in the package directory will run this script


def file_path(string):
    if os.path.isfile(string) or string == "../../" :
        return string
    else:

        raise OSError(string)

def dir_path(string):
    if os.path.isdir(string) or string == "../../" :
        return string
    else:

        raise OSError(string)

def main():
    # create argument parser object
    parser = argparse.ArgumentParser(description = "Quantum Tomography")

    parser.add_argument("-i", "--eval", type = file_path, nargs = 1,
                        metavar = "evalFile", default = None, help = "The full path to the file that contains the data and configuration for the tomography.")

    parser.add_argument("-s", "--save", type = dir_path, nargs = 1,
                        metavar = "outPutFolder", default = None, help = "The full path to the folder where you want the output to be saved. If not included it will not save your data.")

    parser.add_argument("-p", "--pic", action = 'store_true',default = False,
                        help = "Including this will show images of real and imaginary values of the density matrix. If save is also included pictures will be saved and not shown.")


    # parse the arguments from standard input
    args = parser.parse_args()

    try:
        inPutfilePath = args.eval[0]
    except:
        raise ValueError("input not defined")
    try:
        outPutfilePath = args.save[0]
        save = True
    except:
        save = False
    pictures = args.pic

    t = qKLib.Tomography()

    # import the eval file to import both the config and data
    [rho, intensity, fval] = t.importEval(inPutfilePath)
    qKLib.printLastOutput(t)

    if(save and save != "False"):
        if not os.path.exists(outPutfilePath + '/TomoOutPut'):
            os.makedirs(outPutfilePath + '/TomoOutPut')
        # Prints the data to a html file
        FORREPLACE = '<h1 style = "text-align: center;"> Tomography Results</h1>'
        if(pictures):
            import matplotlib.pyplot as plt
            FORREPLACE = FORREPLACE + '<img src = "rhobarReal.png" style = "float: left;" width = "550" height = "auto">'
            FORREPLACE = FORREPLACE + '<img src = "rhobarImag.png" width = "550" height = "auto"><br>'
            qKLib.saveRhoImages(rho,outPutfilePath + '/TomoOutPut')

        FORREPLACE = FORREPLACE + qKLib.matrixToHTML(rho)

        vals = t.getProperties(rho)
        FORREPLACE = str(FORREPLACE) + qKLib.propertiesToHTML(vals)

        # Print out properties of bellSettings
        bs = ''
        if (t.conf['Bellstate'] != 0):
            vals = t.getBellSettings(rho)
            resolution = (np.pi / 2) / (9 * (5 ** 3))
            bs += '<h3>Bell inequality (S_local_realism <= 2)</h3>'
            bs += '<table style = \"width:60%;margin-top:10px;font-size: 15px;padding-bottom:5px;float:none;\"><tr><td style = "font-size: 20;font-weight: 1000;color: rebeccapurple;">Property</td><td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">Value</td>'
            if (vals.shape[1] > 2):
                bs += '<td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">   Error</td></tr>'
            else:
                bs += '<td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;"></td></tr>'
            for v in vals:
                bs += '<tr><td>' + v[0] + '</td><td>' + qKLib.floatToString(v[1]) + ' &deg;</td>'
                if (len(v) > 2):
                    bs += '<td>+/- ' + qKLib.floatToString(v[2]) + '&deg;</td></tr> \n'
                else:
                    bs += "<td></td></tr>"
            bs += '<tr><td>resolution</td><td>' + qKLib.floatToString(resolution * 180 / np.pi) + '&deg;</td><td></tr></tr> \n'
            bs += '</td></tr></table>'
        bs = bs.replace('+/- nan&deg;', 'Error')
        bs = bs.replace('nan &deg;', 'Error')

        # Print out time at the bottom of page
        FORREPLACE = str(FORREPLACE) + '<br><br>\n' + bs
        FORREPLACE += '</font></div>'


        # Update the html file
        fff = "<html><head></head><body>TOREPLACE</body></html>"

        fff = fff.replace('TOREPLACE', str(FORREPLACE))

        with open(outPutfilePath + '/TomoOutPut/outPut.html', 'w') as ff:
            ff.write(fff)
            ff.close()
        if(outPutfilePath[-1] == "\\" or outPutfilePath[-1] == "/"):
            print("Output saved to "+outPutfilePath + 'TomoOutPut')
        else:
            if (outPutfilePath.find("\\") == -1):
                print("Output saved to " + outPutfilePath + '/TomoOutPut')
            else:
                print("Output saved to " + outPutfilePath + '\\TomoOutPut')
    else:
        if(pictures):
            import matplotlib.pyplot as plt
            qKLib.makeRhoImages(rho,plt)
            plt.show()

if __name__ == "__main__":
    main()
