# This script is used to run quantum tomography from the command line.
#
# Running the command Quantum-Tomography from the command line in the package directory will run this script
#
# Arguments:
#
import argparse
import os
import sys                                          #for path to external scripts
import src as qKLib
import numpy as np
from numpy.core.defchararray import add



def file_path(string):
    if os.path.isfile(string) or string == "../../" :
        return string
    else:

        raise FileNotFoundError(string)

def dir_path(string):
    if os.path.isdir(string) or string == "../../" :
        return string
    else:

        raise FileNotFoundError(string)

def main():
    # create argument parser object
    parser = argparse.ArgumentParser(description="Weather Reporter")

    parser.add_argument("-i", "--eval", type=file_path, nargs=1,
                        metavar="evalFile", default= None, help="The full path to the file that contains the data and configuration for the tomography")

    parser.add_argument("-o", "--outPut", type=dir_path, nargs=1,
                        metavar="outPutFolder", default= None, help="The full path to the folder where you want the output to be saved")

    parser.add_argument("-p", "--pic", type=bool, nargs=1,
                        metavar="pictureBool", default=False,
                        help="Setting this to true will provide images of real and imaginary values of the density matrix")
    parser.add_argument("-s", "--save", type=bool, nargs=1,
                        metavar="saveBool", default=False,
                        help=" Setting this to true will save your output data in HTML format")

    # parse the arguments from standard input
    args = parser.parse_args()

    try:
        inPutfilePath = args.eval[0]
    except:
        raise ValueError("input not defined")
    try:
        outPutfilePath = args.outPut[0]
    except:
        raise ValueError("output not defined")
    save = args.save
    pictures = args.pic

    t = qKLib.Tomography()

    # import the eval file to import both the config and data
    [rho, intensity, fval] = t.importEval(inPutfilePath)
    qKLib.printLastOutput(t)

    if(save):
        #Prints the data to a html file
        FORREPLACE = '<h1 style="text-align: center;"> Tomography Results</h1>'
        if(pictures):
            import matplotlib.pyplot as plt
            import warnings
            warnings.filterwarnings('ignore')
            FORREPLACE = FORREPLACE + '<img src="rhobarReal.png" style="float: left;" width="550" height="auto">'
            FORREPLACE = FORREPLACE + '<img src="rhobarImag.png" width="550" height="auto"><br>'
            qKLib.saveRhoImages(rho,outPutfilePath + '/TomoOutPut')
            # if (t.conf['NQubits'] > 1):
            #     fig = plt.figure()
            #     ax1 = fig.add_subplot(111, projection='3d')
            #     xpos = [.5, .5, .5, .5,
            #             1.5, 1.5, 1.5, 1.5,
            #             2.5, 2.5, 2.5, 2.5,
            #             3.5, 3.5, 3.5, 3.5]
            #     ypos = [.5, 1.5, 2.5, 3.5,
            #             .5, 1.5, 2.5, 3.5,
            #             .5, 1.5, 2.5, 3.5,
            #             .5, 1.5, 2.5, 3.5, ]
            #     zpos = np.zeros(16)
            #     dx = [.9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9]
            #     dy = [.9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9]
            #
            #     dz = [.5, 1.5, 2.5, 3.5,
            #           .5, 1.5, 2.5, 3.5,
            #           .5, 1.5, 2.5, 3.5,
            #           .5, 1.5, 2.5, 3.5, ]
            #     dz = rho.flatten().astype(float)
            #
            #     colors = ['r', 'g', 'b', 'purple']
            #
            #     # ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
            #     for i in range(4):
            #         ax1.bar3d(xpos[i * 4:i * 4 + 4], ypos[i * 4:i * 4 + 4], zpos[i * 4:i * 4 + 4], dx[i * 4:i * 4 + 4],
            #                   dy[i * 4:i * 4 + 4], dz[i * 4:i * 4 + 4], edgecolor="black", alpha=.65, color=colors[i])
            #
            #     ax1.axes.set_xticklabels(["|HH>", "|HV>", "|VH>", "|VV>"])
            #     ax1.axes.set_yticklabels(["<HH|", "<HV|", "<VH|", "<VV|"])
            #     ax1.axes.set_xticks([1, 2, 3, 4])
            #     ax1.axes.set_yticks([1, 2, 3, 4])
            #     ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
            #     ax1.axes.set_zlim3d(-1, 1)
            #     plt.title("Rho Real")
            #     plt.savefig(outPutfilePath + 'TomoOutPut/rhobarRreal.png')
            #
            #     fig = plt.figure()
            #     ax1 = fig.add_subplot(111, projection='3d')
            #     xpos = [.5, .5, .5, .5,
            #             1.5, 1.5, 1.5, 1.5,
            #             2.5, 2.5, 2.5, 2.5,
            #             3.5, 3.5, 3.5, 3.5]
            #     ypos = [.5, 1.5, 2.5, 3.5,
            #             .5, 1.5, 2.5, 3.5,
            #             .5, 1.5, 2.5, 3.5,
            #             .5, 1.5, 2.5, 3.5, ]
            #     zpos = np.zeros(16)
            #     dx = [.9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9]
            #     dy = [.9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9
            #         , .9, .9, .9, .9]
            #
            #     dz = [.5, 1.5, 2.5, 3.5,
            #           .5, 1.5, 2.5, 3.5,
            #           .5, 1.5, 2.5, 3.5,
            #           .5, 1.5, 2.5, 3.5, ]
            #     dz = rho.flatten().imag.astype(float)
            #     colors = ['r', 'g', 'b', 'purple']
            #
            #     for i in range(4):
            #         ax1.bar3d(xpos[i * 4:i * 4 + 4], ypos[i * 4:i * 4 + 4], zpos[i * 4:i * 4 + 4], dx[i * 4:i * 4 + 4],
            #                   dy[i * 4:i * 4 + 4], dz[i * 4:i * 4 + 4], edgecolor="black", alpha=.65, color=colors[i])
            #
            #     ax1.axes.set_xticklabels(["|HH>", "|HV>", "|VH>", "|VV>"])
            #     ax1.axes.set_yticklabels(["<HH|", "<HV|", "<VH|", "<VV|"])
            #     ax1.axes.set_xticks([1, 2, 3, 4])
            #     ax1.axes.set_yticks([1, 2, 3, 4])
            #     ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
            #     ax1.axes.set_zlim3d(-1, 1)
            #     plt.title("Rho Imaginary")
            #     plt.savefig(outPutfilePath + 'TomoOutPut/rhobarRimag.png')
            # else:
            #     # 1 Qubit
            #     fig = plt.figure()
            #     ax1 = fig.add_subplot(111, projection='3d')
            #     xpos = [.5, .5,
            #             1.5, 1.5]
            #     ypos = [.5, 1.5,
            #             .5, 1.5]
            #     zpos = np.zeros(4)
            #     dx = [.9, .9,
            #           .9, .9]
            #     dy = [.9, .9,
            #           .9, .9]
            #     dz = [.5, 1.5,
            #           .5, 1.5]
            #     dz = rho.flatten().astype(float)
            #
            #     colors = ['r', 'g', 'b', 'purple']
            #
            #     # ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
            #     for i in range(4):
            #         ax1.bar3d(xpos[i], ypos[i], zpos[i], dx[i], dy[i], dz[i], edgecolor="black", alpha=.65,
            #                   color=colors[i])
            #     ax1.axes.set_xticklabels(["|H>", "|V>"])
            #     ax1.axes.set_yticklabels(["<H|", "<V|"])
            #     ax1.axes.set_xticks([1, 2])
            #     ax1.axes.set_yticks([1, 2])
            #     ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
            #     ax1.axes.set_zlim3d(-1, 1)
            #     plt.title("Rho Real")
            #     plt.savefig(outPutfilePath + 'TomoOutPut/rhobarRreal.png')
            #     # plt.show()
            #
            #     fig = plt.figure()
            #     ax1 = fig.add_subplot(111, projection='3d')
            #     xpos = [.5, .5,
            #             1.5, 1.5]
            #     ypos = [.5, 1.5,
            #             .5, 1.5]
            #     zpos = np.zeros(4)
            #     dx = [.9, .9,
            #           .9, .9]
            #     dy = [.9, .9,
            #           .9, .9]
            #     dz = [.5, 1.5,
            #           .5, 1.5]
            #     dz = rho.flatten().imag.astype(float)
            #
            #     colors = ['r', 'g', 'b', 'purple']
            #
            #     # ax1.bar3d(xpos, ypos, zpos, dx, dy, dz, color='#00ceaa')
            #     for i in range(4):
            #         ax1.bar3d(xpos[i], ypos[i], zpos[i], dx[i], dy[i], dz[i], edgecolor="black", alpha=.65,
            #                   color=colors[i])
            #     ax1.axes.set_xticklabels(["|H>", "|V>"])
            #     ax1.axes.set_yticklabels(["<H|", "<V|"])
            #     ax1.axes.set_xticks([1, 2])
            #     ax1.axes.set_yticks([1, 2])
            #     ax1.axes.set_zticks([-1.0, -.8, -.6, -.4, -.2, 0, .2, .4, .6, .8, 1.0])
            #     ax1.axes.set_zlim3d(-1, 1)
            #     plt.title("Rho Imaginary")
            #     plt.savefig(outPutfilePath + 'TomoOutPut/rhobarRimag.png')

        FORREPLACE = FORREPLACE + qKLib.matrixToHTML(rho)
        # FORREPLACE = FORREPLACE + createtable(rho)

        vals = t.getProperties(rho)
        f = '<h3 >Properties of Rho</h3><table style=\"width:60%;margin-top:10px;font-size: 15px;padding-bottom:5px;float:none;\"><tr><td style="font-size: 20;font-weight: 1000;color: rebeccapurple;">Property</td><td  style="font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">Value</td><td  style="font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">   Error</td></tr>'
        for v in vals:
            f += '<tr><td><div onmouseover="Tip( ' + v[0].replace(" ", "") + 'Tip)" onmouseout="hideTip()">' + v[
                0] + '</td><td name= "' + v[0].replace(" ", "") + '_value">' + str(round(v[1], 5)) + '</td>'
            if (len(v) > 2):
                f += '<td> +/- ' + str(round(v[2], 5)) + '</td></tr>'
            else:
                f += "<td></td></tr>"
        f += "</table>"

        FORREPLACE = str(FORREPLACE) + str(f)

        # Print out properties of bellSettings
        bs = ''
        if (t.conf['Bellstate'] != 0):
            vals = t.getBellSettings(rho)
            resolution = (np.pi / 2) / (9 * (5 ** 3))
            bs += '<h3>Bell inequality (S_local_realism <= 2)</h3>'
            bs += '<table style=\"width:60%;margin-top:10px;font-size: 15px;padding-bottom:5px;float:none;\"><tr><td style="font-size: 20;font-weight: 1000;color: rebeccapurple;">Property</td><td  style="font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">Value</td><td  style="font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">   Error</td></tr>'

            for v in vals:
                bs += '<tr><td>' + v[0] + '</td><td>' + str(round(v[1], 5)) + ' &deg;</td>'
                if (len(v) > 2):
                    bs += '<td>+/- ' + str(round(v[2], 5)) + '&deg;</td></tr> \n'
                else:
                    bs += "<td></td></tr>"
            bs += '<tr><td>resolution</td><td>' + str(resolution * 180 / np.pi) + '&deg;</td><td></tr></tr> \n'
            bs += '</td></tr></table>'
        bs = bs.replace('+/- nan&deg;', 'Error')
        bs = bs.replace('nan &deg;', 'Error')

        # Print out time at the bottom of page
        FORREPLACE = str(FORREPLACE) + '<br><br>\n' + bs
        FORREPLACE += '</font></div>'


        # Update the html file
        fff = "<html><head></head><body>TOREPLACE</body></html>"

        fff = fff.replace('TOREPLACE', str(FORREPLACE))
        if not os.path.exists(outPutfilePath + '/TomoOutPut'):
            os.makedirs(outPutfilePath + '/TomoOutPut')

        with open(outPutfilePath + '/TomoOutPut/outPut.html', 'w') as ff:
            ff.write(fff)
            ff.close()
        print("Output saved to "+outPutfilePath + '\\TomoOutPut')
#
# def createtable(M):
#     s = np.shape(M)
#     res = '<table style=\"border: 1px solid black;border-collapse: collapse;font-size: 15px; table-layout:fixed;width:100%;\">'
#     for i in range(s[0]):
#         res = res+' <tr>'
#         for j in range(s[1]):
#             res = res + '<td style = "border: 1px solid black;">' + str(np.real(M[i,j])) + "<div style=\"color:rebeccapurple;font-weight: bold;display:inline;\">+</div><BR>"+ str(np.imag(M[i,j]))
#             res = res + '<div style=\"color:rebeccapurple;font-weight: bold;display:inline;\">j</div></td>'
#         res = res +'</tr>'
#     res = res+'</table>'
#     d, v = np.linalg.eig(M)
#     eigenVals = "<h5>Eigen Values: "
#     for x in range(0,len(d)):
#         eigenVals = eigenVals+str(round(d[x].real, 5))
#         if(abs(d[x].imag)>.00001):
#             eigenVals = eigenVals+"< div style =\"color:rebeccapurple;font-weight: bold;display:inline;\">+</div>"
#             eigenVals = eigenVals+str(round(d[x].imag, 5))
#             eigenVals = eigenVals+"<div style=\"color:rebeccapurple;font-weight: bold;display:inline;\">j</div>"
#         eigenVals = eigenVals+" , "
#     eigenVals = str(eigenVals)[0:len(str(eigenVals))-2]
#     eigenVals = eigenVals+"</h5>"
#     res = res+eigenVals
#     return res

if __name__ == "__main__":
    main()