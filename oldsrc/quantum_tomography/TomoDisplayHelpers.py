from __future__ import print_function
import numpy as np

"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""


# These are various helper functions used in other main functions in the TomoDisplay file


# This function rounds and converts floats to strings. It handles complex numbers and uses scientific notation
def floatToString(x, html = False):
    if(x == "NA"):
        return x
    if(abs(x.imag) > 10**-8):
        tempIMAG = floatToString(x.imag)
        if(abs(x.real) > 10**-8):
            if(html):
                if(tempIMAG[0] != "-"):
                    return floatToString(x.real) + "<div style = \"color:rebeccapurple;font-weight: bold;display:inline;\"> + i</div>" + tempIMAG
                else:
                    return floatToString(x.real) + "<div style = \"color:rebeccapurple;font-weight: bold;display:inline;\"> - i</div>" + tempIMAG[1:]
            else:
                if (tempIMAG[0] != "-"):
                    return floatToString(x.real) + " + i" + tempIMAG
                else:
                    return floatToString(x.real) + " - i" + tempIMAG[1:]
        else:
            if (html):
                if (tempIMAG[0] != "-"):
                    return "<div style = \"color:rebeccapurple;font-weight: bold;display:inline;\">i</div>" + tempIMAG
                else:
                    return "<div style = \"color:rebeccapurple;font-weight: bold;display:inline;\">-i</div>" + tempIMAG[1:]
            else:
                if (tempIMAG[0] != "-"):
                    return "i" + tempIMAG
                else:
                    return "-i" + tempIMAG[1:]
    else:
        if(x == float("inf")):
            return "inf"
        if isNaN(x):
            return "nan"
        if (abs(x.real) <= 10 ** -8):
            return "0"

        s = "{:e}".format(float(x.real))
        [num, power] = s.split("e")
        if(num[0] == "-"):
            num = num[:5]
        else:
            num = num[:4]

        if(abs(float(power)) > 2):
            return num+"e"+power
        else:
            s = float(num)*10**float(power)
            s = np.around(s, 2-int(power))
            s = str(s)
            return s
# checks if a number is nan
def isNaN(num):
    return num != num
#
# tests = [1*10**-20, 1j*10**-20, 1*10**-20+1j*10**-20, 1*10**-5+1j*10**-20, 1*10**-20+1j*10**-5, 1*10**-5, 1j*10**-5, 0]
# for i in tests:
#     print(str(i)+" : "+ floatToString(i))
