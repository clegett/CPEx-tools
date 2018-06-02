#!/usr/bin/env python3
"""Convert two column spectral data from wavelength to wavenumbers.

Converts spectral data from wavelength to wavenumbers. Data is accepted in
any two column (wavelength | wavenumber)(delimiter)(reflectance, etc...)
format.


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
import itertools
from pint import UnitRegistry
import sys
import csv

__author__ = 'Carey Legett'
__contact__ = 'carey.legett@stonybrook.edu'
__copyright__ = 'Copyright 2018, Carey Legett'
__credits__ = ['Carey Legett']
__date__ = '2018/06/01'
__deprecated__ = False
__email__ = 'carey.legett@stonybrook.edu'
__license__ = 'GPLv3'
__maintainer = 'Carey Legett'
__status__ = 'Development'
__version__ = '0.1a'

ureg = UnitRegistry()


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

    parser = argparse.ArgumentParser(prog='wl_to_wn.py',
                                     description=__doc__,
                                     epilog='''For more information contact 
                                     Chip at carey.legett@stonybrook.edu.''')
    parser.add_argument('-f', '--input-file', required=True,
                        help='''The name of the input file (include path if the
                        file is not in this directory).''', dest='ifile')
    parser.add_argument('-x', '--xunits', type=str, required=True,
                        help='''The units that the x data are in (e.g. 
                        'micrometers' or 'wavenumbers'.''', dest='xunits')
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument('-d', '--delimiter', help='''The delimiter for the input
                       file. A single character representing the actual 
                       delimiter. Mutually exclusive with -w or 
                       --whitespace''', dest='in_delim')
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
    parser.add_argument('-s', '--output-delimiter', required=False, help='''The
                        delimiter to be used in the output file. Defults to tab
                        .''', dest='out_delim')
    parser.add_argument('-v', '--verbose', action='store_true', dest='v',
                        help='Give status information during processing.')

    return parser.parse_args()


def read_input_file(filename, delimiter, whitespace_flag, skip_lines, xunits,
                    verbose):
    if verbose:
        print('Reading input file...')

    try:
        with open(filename, 'rt', newline='') as infile:
            if whitespace_flag:
                rawdata = [line.split() for line in
                           itertools.islice(infile, skip_lines, None)]
            else:
                rawdata = [line.split(delimiter) for line in
                           itertools.islice(infile, skip_lines, None)]
        if not rawdata:
            sys.exit('No data in file {}'.format(filename))
    except IOError as e:
        sys.exit('I/O error: file {}: {}'.format(filename, e))

    mydata = [[float(astring) for astring in inner] for inner in rawdata]
    print(mydata)
    if verbose:
        print('Assigning units.')
    mydata = [[row[0] * ureg[f'{xunits}'], row[1]] for row in mydata]
    # If this fails, it is likely that the data read from the input file is not
    #  numeric. It could also happen if the delimiter was not correctly
    # specified.

    if verbose:
        print('Done with input.')

    return mydata


def write_to_output(somedata, out_filename, out_delim, verbose_flag):
    if verbose_flag:
        print('Writing output file')
        print(somedata)
    try:
        with open(out_filename, 'w', newline='') as outfile:
            writer = csv.writer(outfile, delimiter=out_delim)
            header = ['wn', 'data']

            writer.writerow(header)
            writer.writerows(somedata)
    except IOError as e:
        sys.exit(f'Output file error: {e}, {sys.exc_info()}')


if __name__ == '__main__':
    args = parse_cmdln_args()
    x = ureg[f'{args.xunits}']
    data = read_input_file(args.ifile, args.in_delim, args.whitespace,
                           args.skip, args.xunits, args.v)

    #print(data)
    data = [[line[0].to(1 / ureg.centimeter, 'sp').magnitude, line[1]] for line in data]

    write_to_output(data, args.ofile, args.out_delim, args.v)
