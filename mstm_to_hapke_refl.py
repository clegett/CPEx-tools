#!/usr/bin/env python3

#
# This file reads MSTM output files and outputs a CSV after calculating
#  the reflectance from the MSTM data.
#
# MSTM files should be named [wavelength].dat (e.g 700nm.dat)
# All input files to be read should be in the directory this
#  file is run from. This program will ignore commandline args.
#
# This program was written by Chip Legett carey.legett@stonybrook.edu
# 17 November 2015
#
# Updates:
#  18 Mar 2017: added B&H phase function and detection and handling of 
#               fixed vs random orientation runs. -CL
#


import sys   # for sys.exit file read failure
import glob  # for getting list of files in directory
import re    # regular expressions, because fml
import os    # path handling for filenames
import csv   # csv handling for output file
import math  # we want access to pi and cosine
from operator import itemgetter  # used to sort the list for output


# Later on, we're going to need an H function, it is defined here
def h(mu, gamma):
    top = 1 + (2 * mu)
    bottom = 1 + (2*gamma*mu)
    return top/bottom


# I'm going to have the reflectance calculation be a function just to
#  keep it a little cleaner down below in the loop where it is used
def refl(ssa, ideg, edeg, phase=1.0, backscatter=0.0):
    """Calculate the Hapke Bidirectional Reflectance Function

    Keyword arguments:
    ssa -- the single scattering albedo
    ideg -- the incidence angle in degrees
    edeg -- the emergence angle in degrees
    phase -- the phase function (default 1.0)
    backscatter -- the backscatter function (default 0.0)
    """

    mu0 = math.cos(math.radians(ideg))
    mu = math.cos(math.radians(edeg))
    gamma = math.sqrt(1.0 - ssa)
    hmu0 = h(mu0, gamma)
    hmu = h(mu, gamma)
    refl = (ssa / (4 * math.pi)) * (mu0 / (mu0 + mu)) * (
            (1 + backscatter) * phase + (hmu0 * hmu) - 1)
    return refl


def run_script(a_dir='./'):
    # get list of .dat files in directory
    file_list = glob.glob('{}*nm.dat'.format(a_dir))

    # check and see if there is anything in the file list and exit if not
    if len(file_list) == 0:
        sys.exit('No input files found in directory {}!'.format(a_dir))

    #
    # The data list is going to be formatted like this:
    #
    # index  0       1    2    3    4    5   6    7    8
    # 0[[wavelength,Qext,Qabs,Qsct,Asym,SSA,S11,Refl,flag],
    # 1 [wavelength,Qext,Qabs,Qsct,Asym,SSA,S11,Refl,flag],
    # n...
    # ]
    #
    # where data[1][3] is Qsct for the second row
    #
    # Note that this is a list of lists.
    #
    # We'll declare the variable here, with one row for every file
    #  and 9 columns for the 9 data points extracted from each file
    #

    data = [[0 for i in range(9)] for j in range(len(file_list))]

    # I'm also going to add a counter variable to iterate over
    #  to help me keep track of which row of the array I'm writing to

    row = 0

    #
    # Now let's get the files open and collect some data
    #

    # I'll declare the s11index and ext_search_string variables to the fixed
    # orientation values primarily for scope here. These should be detected
    # and set in the loop below before they are used.

    s11index = 2
    ext_search_string = 'unpolarized total ext'

    # In order to calculate the phase function we need to know the geometric
    #  cross section (G) of the particle.

    # This blurb of code makes the input function usable across Python versions
    try:
        input = raw_input
    except NameError:
        pass

    print('\nThis will be used to calculate the geometric cross section:')
    print('What is the RADIUS of the enclosing sphere?')
    user_input = input()
    enc_radius = float(user_input.strip(' \t\n\r'))
    G = math.pi * enc_radius**2

    # Loop on all .dat files in directory:
    for currentfile in file_list:
        # open file
        with open(currentfile) as f:
            # get wavelength from filename
            data[row][0] = float(os.path.splitext(f.name)[0][0:-2])
            # This will be our super simple state machine switch variable
            state = 0
            # search through file until we find the line before the
            # efficiencies
            for line in f:
                # remove leading and trailing whitespace
                clean_line = line.strip()
                # State == 1 means we have determined that this line contains
                #  efficiency data and we're going to read it
                if state == 1:
                    state = 0
                    # split string on spaces
                    elements = clean_line.split()
                    # get Qext
                    data[row][1] = float(elements[0])
                    # get Qabs
                    data[row][2] = float(elements[1])
                    # get Qsct
                    data[row][3] = float(elements[2])
                    # get Asym
                    data[row][4] = float(elements[3])
                    # calculate SSA
                    data[row][5] = data[row][3] / data[row][1]
                    # If the SSA is > 1, we change it to 1 so we don't
                    #  break a sqrt later on, if we do this, we set the
                    #  "flag" parameter to 1 so we can make note of it
                    #  in the output file later on
                    if data[row][5] > 1:
                        data[row][5] = 1
                        data[row][8] = 1
                    continue
                # State == 2 means that the next line contains the 150 degree
                #  scattering matrix data and we're going to read it
                elif state == 2:
                    # split string on spaces
                    matrix = clean_line.split()
                    # get S11
                    data[row][6] = float(matrix[s11index])
                    # Calculate phase function (See Bohren and Huffman pp
                    # 383-384)
                    k = (2 * math.pi)/data[row][0]
                    dcsca = data[row][6] / k**2
                    Csca = data[row][3] * G
                    phase = 4*math.pi*(1/Csca)*dcsca
                    # We now have everything necessary to calculate the
                    # reflectance
                    data[row][7] = refl(data[row][5], 30, 0, phase)
                    # reset state machine
                    state = 0
                    # increment row for next file
                    row += 1
                    # break out of this file
                    break
                # Here we want to search the line for the text "random
                # orientation calculations" to know what type of model run
                # this is. This is necessary to determine the correct format of
                # the scattering matrix data used a little further down in
                # this loop
                elif (re.search('random orientation calculations', clean_line)
                      is not None):
                    # This file was output from a random orientation run
                    #  This means that the scattering matrix data is formatted:
                    #   theta 11 12 13 14 22 23 24 33 34 44
                    #  we want 11, so we'll grab the second field in a minute
                    s11index = 1
                    ext_search_string = "total ext"
                    continue
                # Here we want to search the line for the text "fixed
                # orientation calculations" to know what type of model run
                # this is. This is necessary to determine the correct format
                # of the scattering matrix data used a little further down
                # in this loop
                elif (re.search('fixed orientation calculations', clean_line)
                      is not None):
                    # This file was output from a fixed orientation run
                    # This means that the scattering matrix data is formatted:
                    # theta phi 11 12 13 14 21 22 23 24 31 32 33 34 41 42 43 44
                    # we want 11, so we'll grab the third field in a minute
                    s11index = 2
                    ext_search_string = 'unpolarized total ext'
                    continue
                # Here we want to search this line for the phrase "total ext"
                #  if the line contains this text, the next line contains the
                #  efficiency data and we want to catch that on the next loop
                elif re.search(ext_search_string, clean_line) is not None:
                    # The next line contains efficiency data, get ready to
                    # catch it
                    state = 1
                    # start the next iteration of the for loop
                    continue
                # Here we want to search the line for "149.00" because this
                # means that the next line contains the scattering matrix
                # information for a scattering angle of 150 degrees and we want
                # to catch that on the next loop
                elif re.search('149.00', clean_line) is not None:
                    # The next line contains 150 deg scattering matrix
                    state = 2
                    # start the next iteration of the for loop
                    continue

    # sort the data list by wavelength
    data.sort(key=itemgetter(0))
    # write array to output file
    return data


def write_data(data, some_dir='./'):
    with open('{}output.csv'.format(some_dir), 'w') as out:
        writer = csv.writer(out)
        writer.writerow(['wavelength', 'Qext', 'Qabs', 'Qsct', 'Asym', 'SSA',
                         'S11', 'Refl', 'flag'])
        for values in data:
            writer.writerow(values)


if __name__ == '__main__':
    a_directory = './'
    my_data = run_script(a_directory)
    write_data(my_data, a_directory)
