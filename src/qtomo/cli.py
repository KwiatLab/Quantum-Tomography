import argparse
import os

import numpy as np

import qtomo
from qtomo.tomography import Tomography


def main() -> None:
    parser = argparse.ArgumentParser(description="Quantum Tomography")

    _ = parser.add_argument(
        "filename"
        type=str,
        nargs=1,
        metavar="FILE",
        required=True,
        help="Path to the tomography data file.",
    )

    _ = parser.add_argument(
        "-c",
        "--config",
        type=str,
        nargs=1,
        metavar="CONFIG",
        default=None,
        help="Path to the config file to use. If not specified, the default configuration will be used.",
    )

    _ = parser.add_argument(
        "-w",
        "--write",
        type=str,
        nargs=1,
        metavar="OUTDIR",
        default=None,
        help="Path directory to save tomography data. If not specified, data will be displayed, but not saved.",
    )

    _ = parser.add_argument(
        "--skip-real-imag",
        action="store_true",
        default=False,
        help="Whether to incldue the skip showing the real and imaginary density matrix components and just show the overall plot.",
    )

    args = parser.parse_args()

    args



    t = Tomography()

    # import the eval file to import both the config and data

    if save and save != "False":
        if not os.path.exists(outPutfilePath + "/TomoOutPut"):
            os.makedirs(outPutfilePath + "/TomoOutPut")
        # Prints the data to a html file
        FORREPLACE = '<h1 style = "text-align: center;"> Tomography Results</h1>'
        if pictures:
            import matplotlib.pyplot as plt

            FORREPLACE = (
                FORREPLACE + '<img src = "rhobarReal.png" style = "float: left;" width = "550" height = "auto">'
            )
            FORREPLACE = FORREPLACE + '<img src = "rhobarImag.png" width = "550" height = "auto"><br>'
            qtomo.saveRhoImages(rho, outPutfilePath + "/TomoOutPut")

        FORREPLACE = FORREPLACE + qtomo.matrixToHTML(rho)

        vals = t.getProperties(rho)
        FORREPLACE = str(FORREPLACE) + qtomo.propertiesToHTML(vals)

        # Print out properties of bellSettings
        bs = ""
        if t.conf["Bellstate"] != 0:
            vals = t.getBellSettings(rho)
            resolution = (np.pi / 2) / (9 * (5**3))
            bs += "<h3>Bell inequality (S_local_realism <= 2)</h3>"
            bs += '<table style = "width:60%;margin-top:10px;font-size: 15px;padding-bottom:5px;float:none;"><tr><td style = "font-size: 20;font-weight: 1000;color: rebeccapurple;">Property</td><td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">Value</td>'
            if vals.shape[1] > 2:
                bs += '<td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;">   Error</td></tr>'
            else:
                bs += (
                    '<td  style = "font-size: 20;font-weight: 1000;color: rebeccapurple;padding-bottom:5px;"></td></tr>'
                )
            for v in vals:
                bs += "<tr><td>" + v[0] + "</td><td>" + qtomo.floatToString(v[1]) + " &deg;</td>"
                if len(v) > 2:
                    bs += "<td>+/- " + qtomo.floatToString(v[2]) + "&deg;</td></tr> \n"
                else:
                    bs += "<td></td></tr>"
            bs += (
                "<tr><td>resolution</td><td>"
                + qtomo.floatToString(resolution * 180 / np.pi)
                + "&deg;</td><td></tr></tr> \n"
            )
            bs += "</td></tr></table>"
        bs = bs.replace("+/- nan&deg;", "Error")
        bs = bs.replace("nan &deg;", "Error")

        # Print out time at the bottom of page
        FORREPLACE = str(FORREPLACE) + "<br><br>\n" + bs
        FORREPLACE += "</font></div>"

        # Update the html file
        fff = "<html><head></head><body>TOREPLACE</body></html>"

        fff = fff.replace("TOREPLACE", str(FORREPLACE))

        with open(outPutfilePath + "/TomoOutPut/outPut.html", "w") as ff:
            ff.write(fff)
            ff.close()
        if outPutfilePath[-1] == "\\" or outPutfilePath[-1] == "/":
            print("Output saved to " + outPutfilePath + "TomoOutPut")
        else:
            if outPutfilePath.find("\\") == -1:
                print("Output saved to " + outPutfilePath + "/TomoOutPut")
            else:
                print("Output saved to " + outPutfilePath + "\\TomoOutPut")
    else:
        if pictures:
            import matplotlib.pyplot as plt

            qtomo.makeRhoImages(rho, plt)
            plt.show()
