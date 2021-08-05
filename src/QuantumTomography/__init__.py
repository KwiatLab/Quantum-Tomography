import sys
from .TomoClass import *
from .TomoDisplay import *
"""
Copyright 2020 University of Illinois Board of Trustees.
Licensed under the terms of an MIT license
"""


"""CHECK OUT THE REFERENCE PAGE ON OUR WEBSITE :
https://quantumtomo.web.illinois.edu/Doc/"""



'''
    Comments Should be formatted like the following:
    (this comment is written in single quotes so that it is not picked up by the documentation generator).

    function(parameter 1, parameter 2)
    Desc: This is where the description is written.
	
    Parameters
    ----------      (10 hyphens, one for each letter)
    parameter 1 : type
    description
    parameter 2 : type
    description
                    (make sure the white spaces between are deleted all the way to the left side of the screen)
    Returns
    -------     (7 hyphens, one for each letter)
    returnValue 1 : type
    description
    returnValue 2 : type
    description
	2
    See Also
     ------   (space, 6 hyphens, another space)
    function1;function2;function2      (these should be the functionTitle without parentheses, as it appears on the Table of Contents
'''


# # # # # # # # # #  #
# # TO THE DEVELOPER #
# # # # # # # # # #  #

# Everything in this code is based on the Matlab files you can download on our tomography website, a more detailed
# comments and definations of variables could be found in those m files, The function names are as same as possible.
#
# 'UseDerivatives' is used automatically when doing 1 detctor per qubit, and is not used when doing 2 detectors per qubit.
# This cause the program to run faster and provide more accurate results in both cases
