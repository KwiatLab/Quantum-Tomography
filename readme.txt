
######################
## TO THE DEVELOPER ##
######################

# Everything in this code is based on the Matlab files you can download on our tomography website, a more detailed
# comments and definations of variables could be found in those m files, The function names are as same as possible.

# If you are debugging this program there is a section underneath the sample sets that involve some usefull debugging
# data sets. You can set the state or choose to have it be random and it will create data.
#
# 'UseDerivatives' is used automatically when doing 1 detctor per qubit, and is not used when doing 2 detectors per qubit.
# This cause the program to run faster and provide more accurate results in both cases
#
#
#   Drift Correction is just set off by default. and there is no implementation of it in this code.
#   I dont know the equation to properly implement it