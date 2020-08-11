#!/Users/maurik/P3/bin/python3
#
# This is an example script showing you how you can create a MiniDST file that has a precisely selected
# set of variables in the output.
#
import sys
import os
import re
from ctypes.util import find_library
import ROOT as R

if 'DYLD_LIBRARY_PATH' in os.environ:
#     print(os.environ['DYLD_LIBRARY_PATH'])
    pass
else:
     print("'DYLD_LIBRARY_PATH' is not set.")

root_lib_fail = R.gSystem.Load("libLcioReader")  # Use ROOT to load the LcioReader
if root_lib_fail:
     print("I could not load the libLcioReader into Python using ROOT.\n")
     print("This cab be caused by several different issues: \n")
     print("1 - Your LD_LIBRARY_PATH (DYLD_LIBRARY_PATH on MacOS) is not set to the directory with the libraries\n")
     print("2 - Your operating system strips the LD_LIBRARY_PATH from the environment of Python when it starts up.\n")
     print("    This will give a warning above. For case2 try changing the first line of this code to use the\n")
     print("    correct version of python, or start the script with 'python script.py arguments\n")
     print("")
     sys.exit()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser("Selecting MiniDST writer. ")
    parser.add_argument("--debug", "-d", action="count", help="Increase debug level.", default=0)
    parser.add_argument("--nevent", "-n", type=int, help="Set the number of events to process.", default=0)
    parser.add_argument("--outfile", "-o", type=str, help="Set the output file name", default="minidst.root")
    parser.add_argument('inputfiles', metavar='infile', type=str, nargs='+', help='input files')

    args = parser.parse_args(sys.argv[1:])
    # print(args.inputfiles)

    reader = R.LcioReader(args.inputfiles)
    reader.SetDebugLevel(args.debug)
    reader.SetOutputFileName(args.outfile);
    reader.DefineBranchMap()  # This sets up the names for the branches, but does not add them to the TTree
    column_list = [x.c_str() for x in reader.GetBranchNames()]
    if args.debug > 1:
        print("Column list:")
        print(column_list)

    #
    # Here you can add your own selection of what branches to deactivate.
    #
    deactivate = [ x for x in column_list if re.match("v0_", x) or re.match("svt_", x)]  # Filter out all v0_ and svt_
    deactivate.append("track_covmatrix")

    for name in deactivate:
        reader.SetBranchActive(name, False)  # turn off all the selected branches.

    if args.debug:
        print("Only writing out:")
        column_list = [x.c_str() for x in reader.GetActiveBranchNames()]
        print(column_list)

    reader.SetBranchMap();
    reader.Run(args.nevent);
    reader.End();

