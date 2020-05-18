# This script is used to run quantum tomography from the command line.
#
# Running the command Quantum-Tomography from the command line in the package directory will run this script
#
# Arguments:
#
import argparse
# from .utils import *
import os

def dir_path(string):
    if os.path.isfile(string):
        return string
    else:
        raise FileNotFoundError(string)

def main():
    # create argument parser object
    parser = argparse.ArgumentParser(description="Weather Reporter")

    parser.add_argument("-q", "--eval", type=dir_path, nargs=1,
                        metavar="evalFile", default= None, help="the full path to the file that contains the data and configuration for the tomography")


    # parse the arguments from standard input
    args = parser.parse_args()

    print(args.eval[0])
    print(str(type(args.eval[0])))

if __name__ == "__main__":
    main()