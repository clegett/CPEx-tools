#!/usr/bin/env python3
"""Removes the continuum (background) of data over a specified range.

Removes the continuum of specified x,y data from a delimited text file by 
calculating a linear fit to two points provided by the user and normalizing 
the data such that the linear fit has a value of 1. Output writes to a CSV file
with two to four columns: "x, [y,] [continuum,] continuum_removed_y" as
controlled by command line arguments.

This program is free software: you can redistribute it and/or modify it under
the terms of the GNU General Public License as published by the Free Software
Foundation, either version 3 of the License, or (at your option) any later
version.

This program is distributed in the hope that it will be useful, but WITHOUT
ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or FITNESS
FOR A PARTICULAR PURPOSE. See the GNU General Public License for more details.

You should have received a copy of the GNU General Public License along with
this program. If not, see <http://www.gnu.org/licenses/>.
"""

import argparse
import sys
import itertools
import csv
import numpy as np

__author__ = 'Carey Legett'
__contact__ = 'carey.legett@stonybrook.edu'
__copyright__ = 'Copyright 2018, Stony Brook University'
__credits__ = ['Carey Legett', 'Freenode #r/linux', 'Freenode #python']
__date__ = '2018/02/19'
__deprecated__ = False
__email__ = 'carey.legett@stonybrook.edu'
__license__ = 'GPLv3'
__maintainer = 'Carey Legett'
__status__ = 'Development'
__version__ = '1.0'


def parse_cmdln_args():
    """A function to parse command line arguments for this module.
    
    This function is used when this module is run as a script to parse
    commandline arguments into the appropriate variables. It also provides
    useage help information when run with the -h flag.
    
    Example:
        To invoke the help text and get information about the proper command
        line flags, run the script with the -h or --help flags.

        $ python3 continuum_removal.py -h
        or
        $ python3 continuum_removal.py --help
    
    Returns:
        args (argparse.Namespace): An object containing the state of the
        variables set via command line options. Individual elements are
        accessible with standard `args.element`_ notation.
    
    """

    parser = argparse.ArgumentParser(prog='continuum_removal.py',
                                     description=__doc__,
                                     epilog='''For more information contact 
                                     Chip at carey.legett@stonybrook.edu.''')
    parser.add_argument('-f', '--input-file', required=True,
                        help='''The name of the input file (include path if the
                        file is not in this directory).''', dest='ifile')
    parser.add_argument('-a', '--xstart', type=int, required=True,
                        help='''The x value at which to begin the continuum 
                        calculation.''', dest='xa')
    parser.add_argument('-b', '--xstop', type=int, required=True, help='''The x
                        value at which to end the continuum calculation.''',
                        dest='xb')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--delimiter', help='''The delimiter for the input
                       file. A single character representing the actual 
                       delimiter. Mutually exclusive with -w or 
                       --whitespace''', dest='delim')
    group.add_argument('-w', '--whitespace', help='''Indicates the delimiter is
                       either tabs or spaces. Mutually exclusive with -d or 
                       --delimiter.''', action='store_true', dest='whitespace')
    parser.add_argument('-n', '--skip-lines', default=0, help='''Skip this many
                        lines of the input file before reading in data.''',
                        type=int, dest='skip')
    parser.add_argument('-o', '--output-file', nargs='?', default='output.txt', 
                        help='''The name of the output file (include path 
                        if you want the file somewhere other than in the 
                        current directory). Defaults to output.txt.''',
                        dest='ofile')
    parser.add_argument('-i', '--include-input', action='store_true', dest='i',
                        help='''Write the original input data before the
                        (optional) continuum fit and continuum removed 
                        data.''')
    parser.add_argument('-c', '--continuum-output', action='store_true',
                        dest='c', help='''Write the continuum fit for each 
                        input x value after the (optional) original input data 
                        and before the continuum removed data.''')
    parser.add_argument('-v', '--verbose', action='store_true', dest='v', 
                        help='Give status information during processing.')
    
    return parser.parse_args()


def read_input_file(input_filename, input_delimiter, whitespace_flag, 
                    lines_to_skip, verbose_flag):
    """Reads the specified input file and returns a list of the data therein.

    This function reads an input file using a specified delimiter, skipping a
    given number of lines before reading. It also accepts a verbosity flag to
    enable some progress text while reading.

    Args:
        input_filename (str): The file to be read
        input_delimiter (str): The delimiter used in the input file. Set this
            to `None` if the delimiter is some form of whitespace. 
        whitespace_flag (bool): True if the delimiter is some form of
            whitespace, False otherwise. Must be false if `input_delimiter` is
            not `None`.
        lines_to_skip (int): The number of lines to skip at the beginning of
            the file before ingesting data. 
        verbose_flag (bool): True for progress output, False otherwise.

    Returns:
        data ([float,float]): A list of lists of floats containing the data
            read in from the input file.

    Raises:
        IOError: if there is a problem opening or reading the input file
        ValueError: if the data that is read from the file cannot be converted
            to a float. (e.g. the read data is not numeric)

    """
    if verbose_flag:
        print('Reading input file...')

    # Check for the version of python being used and use appropriate flags for
    # opening the input file as necessary
    if sys.version_info[0] == 2:
        raccess = 'rb'
        kwargs = {}
    else:
        raccess = 'rt'
        kwargs = {'newline': ''}
    
    try:
        with open(input_filename, raccess, **kwargs) as infile:
            if whitespace_flag:
                rawdata = [line.split() for line in itertools.islice(infile,
                           lines_to_skip, None)]
            else:
                rawdata = [line.split(input_delimiter) for line in
                           itertools.islice(infile, lines_to_skip, None)]
        if not rawdata:
            sys.exit('No data in file {}'.format(input_filename))
    except IOError as e:
        sys.exit('I/O error: file {}: {}'.format(input_filename, e))

    mydata = [[float(astring) for astring in inner] for inner in rawdata]
    # If this fails, it is likely that the data read from the input file is not
    #  numeric. It could also happen if the delimiter was not correctly
    # specified.

    return mydata


def remove_continuum(mydata, xa, xb, verbose_flag):
    """A function to remove the linear continuum between two points.
    
    This function takes a list of lists of floats containing x,y data, a
    starting x value, and an ending x value and calculates a linear fit to the
    data at those two points, extrapolated over the provided domain of x. The
    data is then normalized such that the calculated continuum equals 1. This
    can be used to prepare data for analyses such as measuring the band depth
    of absorptions in spectroscopic data.
    
    Args:
        mydata ([float, float]): A list of lists of two floats representing the
            unnormalized initial data.
        xa (float): The value of x at which to start the linear fit. This
            must be an x value in the range of the data list since the
            searching doesn't check to make sure xa falls within the bounds of
            the data.
        xb (float): The value of x at which to stop the linear fit. This
            must be an x value in the range of the data list since the
            searching doesn't check to make sure xa falls within the bounds of
            the data.
        verbose_flag (bool): True for progress output, False otherwise.
        
    Returns:
        output ([float, float]): A list of lists of two to four floats. The
            first element in each inner list will be the input x value. The 
            second element will be the input y value.  The third element will 
            be the continuum value calculated for that x. The last element 
            will be the continuum-removed y data.
    
    """
    if verbose_flag:
        print('Calculating continuum line')

    start_index = np.abs(np.array(mydata)[:, 0] - xa).argmin()
    stop_index = np.abs(np.array(mydata)[:, 0] - xb).argmin()

    startx = mydata[start_index][0]
    stopx = mydata[stop_index][0]

    starty = mydata[start_index][1]
    stopy = mydata[stop_index][1]
    
    if verbose_flag: 
        print('Start: ({!s},{!s})'.format(startx, starty))
        print('Stop: ({!s},{!s})'.format(stopx, stopy))
    
    slope = (stopy-starty)/(stopx-startx)
    intercept = starty-(slope*startx)
    
    if verbose_flag:
        print('Continuum equation: y={}*x+{}'.format(slope, intercept))
        print('Removing continuum')
    
    this_output = []
    for element in mydata:
        continuum = slope*element[0]+intercept
        crdata = element[1]/continuum
        temp_output = [element[0], element[1], continuum, crdata]
        this_output.append(temp_output)

    return this_output


def write_cont_removed_output(myoutput, output_filename, write_input_flag,
                              write_cont_flag, verbose_flag):
    """A function to write the provided list of lists of floats to a CSV file.
    
    This function writes the continuum removed data to a CSV file. Depending on
    the arguments passed, it may also write the original input data, and the
    calculated continuum for each x value.
    
    Args:
        myoutput ([float, float, float, float]): A list of lists of four floats
            that will be written, in part or in whole to the output file.
        output_filename (str): The name of the output file.
        write_input_flag (bool): True if we want to include the original data
            in the output file, False otherwise.
        write_cont_flag (bool): True if we want to include the value of the
            continuum at each x value in the output file, False otherwise.
        verbose_flag (bool): True for progress message, False otherwise.
    
    Raises:
        IOError: For any problem opening or writing to the output file.
    
    Output:
        A CSV file containing two to four columns of data as selected with the
            arguments passed to the function.

"""
    if verbose_flag:
        print('Writing output file')

    # Check for the version of python being used and use appropriate flags for
    # opening the input file as necessary
    if sys.version_info[0] == 2:
        waccess = 'wb'
        kwargs = {}
    else:
        waccess = 'w'
        kwargs = {'newline': ''}

    try:
        with open(output_filename, waccess, **kwargs) as outfile:
            writer = csv.writer(outfile)
            header = ['x']
            filtered_output = []

            for line in myoutput:
                temp = [line[0]]
                if write_input_flag:
                    temp.append(line[1])
                if write_cont_flag:
                    temp.append(line[2])
                temp.append(line[3])
                filtered_output.append(temp)

            if write_input_flag:
                header.append('original')
            if write_cont_flag:
                header.append('continuum')
            header.append('cont_removed')

            writer.writerow(header)
            writer.writerows(filtered_output)
    except IOError as e:
        sys.exit('Output file error: {}, {}'.format(e, sys.exc_info()))


if __name__ == '__main__':
    args = parse_cmdln_args()
    data = read_input_file(args.ifile, args.delim, args.whitespace, args.skip,
                           args.v)
    output = remove_continuum(data, float(args.xa), float(args.xb), args.v)
    write_cont_removed_output(output, args.ofile, args.i, args.c, args.v)
