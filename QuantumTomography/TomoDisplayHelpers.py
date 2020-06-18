"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""

__author__ = 'Quoleon/Turro'
"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE : http://research.physics.illinois.edu/QI/Photonics/Quantum-Tomography_lib_Ref/"""

import numpy as np


#These are various helper functions used in other main functions in the TomoDisplay file


#This function rounds and converts floats to strings. It handles complex numbers and uses scientific notation
def floatToString(x,html=False):
    if(x == "NA"):
        return x
    if(abs(x.imag) > 10**-8):
        if(abs(x.real) > 10**-8):
            if(html):
                return floatToString(x.real) + "<div style=\"color:rebeccapurple;font-weight: bold;display:inline;\">+<BR>j</div>" + floatToString(x.imag)
            return floatToString(x.real) + " + i" + floatToString(x.imag)
        else:
            if (html):
                return "<div style=\"color:rebeccapurple;font-weight: bold;display:inline;\">+<BR>j</div>" + floatToString(x.imag)
            return "i" + floatToString(x.imag)
    else:
        if(x==float("inf")):
            return "inf"
        if isNaN(x):
            return "nan"
        if (abs(x.real) <= 10 ** -8):
            return "0"

        s = "{:e}".format(float(x.real))
        [num,power] = s.split("e")
        num = num[:4]

        if(abs(float(power)) > 2):
            return s
        else:
            s = float(num)*10**float(power)
            return str(s)
#checks if a number is nan
def isNaN(num):
    return num != num
#
# tests = [1*10**-20,1j*10**-20,1*10**-20+1j*10**-20,1*10**-5+1j*10**-20,1*10**-20+1j*10**-5,1*10**-5,1j*10**-5,0]
# for i in tests:
#     print(str(i)+" : "+ floatToString(i))