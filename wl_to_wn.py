#!/usr/bin/env python3
"""Convert two column spectral data from wavelength to wavenumbers.

Converts spectral data from wavelength to wavenumbers. Data is accepted in
any two column (wavelength | wavenumber)(delimiter)(reflectance, etc...) format.


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
from enum import Enum
from pint import UnitRegistry

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

HC = 1.23984193


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
    parser.add_argument('-o', '--output-file', nargs='?', default='output.txt',
                        help='''The name of the output file (include path 
                        if you want the file somewhere other than in the 
                        current directory). Defaults to output.txt.''',
                        dest='ofile')
    parser.add_argument('-v', '--verbose', action='store_true', dest='v',
                        help='Give status information during processing.')

    return parser.parse_args()


class SpectralUnits(Enum):
    WAVENUMBER = 0
    MICRON = 1
    NANOMETER = 2
    ANGSTROM = 3
    EV = 4
    KEV = 5


class SpectralData(object):
    def __init__(self, xdata=None, ydata=None, xunits=None, yunits=None,
                 name=None):
        if xdata is None:
            self.xdata = []
        else:
            self.xdata = xdata
        if ydata is None:
            self.ydata = []
        else:
            self.ydata = ydata
        self.xunits = xunits
        self.yunits = yunits
        self.name = name

    def convert_xdata(self, to_units, from_units=None):
        if from_units is None:
            from_units = self.xunits
        if from_units is to_units:
            pass
        elif from_units is SpectralUnits.WAVENUMBER:
            if to_units is SpectralUnits.MICRON:
                self.xdata = [1e4 / x for x in self.xdata]
                self.xunits = to_units
            elif to_units is SpectralUnits.NANOMETER:
                self.xdata = [1e7 / x for x in self.xdata]
                self.xunits = to_units
            elif to_units is SpectralUnits.ANGSTROM:
                self.xdata = [1e8 / x for x in self.xdata]
                self.xunits = to_units
            elif to_units is SpectralUnits.EV:
                self.xdata = [HC / 1e-4 * x for x in self.xdata]
                self.xunits = SpectralUnits.EV
            elif to_units is SpectralUnits.KEV:
                self.xdata = [HC / 1e-7 * x for x in self.xdata]
                self.xunits = SpectralUnits.KEV
            else:
                pass
        elif from_units is SpectralUnits.MICRON:
            if to_units is SpectralUnits.WAVENUMBER:
                self.xdata = [1e4 / x for x in self.xdata]
                self.xunits = SpectralUnits.WAVENUMBER
            elif to_units is SpectralUnits.NANOMETER:
                self.xdata = [1e3 * x for x in self.xdata]
                self.xunits = SpectralUnits.NANOMETER
            elif to_units is SpectralUnits.ANGSTROM:
                self.xdata = [1e4 * x for x in self.xdata]
                self.xunits = SpectralUnits.ANGSTROM
            elif to_units is SpectralUnits.EV:
                self.xdata = [HC / x for x in self.xdata]
                self.xunits = SpectralUnits.EV
            elif to_units is SpectralUnits.KEV:
                self.xdata = [HC / x / 1e-3 for x in self.xdata]
                self.xunits = SpectralUnits.KEV
            else:
                pass
        elif from_units is SpectralUnits.NANOMETER:
            if to_units is SpectralUnits.WAVENUMBER:
                self.xdata = [1e7 / x for x in self.xdata]
                self.xunits = SpectralUnits.WAVENUMBER
            elif to_units is SpectralUnits.MICRON:
                self.xdata = [x / 1e3 for x in self.xdata]
                self.xunits = SpectralUnits.MICRON
            elif to_units is SpectralUnits.ANGSTROM:
                self.xdata = [1e1 * x for x in self.xdata]
                self.xunits = SpectralUnits.ANGSTROM
            elif to_units is SpectralUnits.EV:
                self.xdata = [HC * 1e3 / x for x in self.xdata]
                self.xunits = SpectralUnits.EV
            elif to_units is SpectralUnits.KEV:
                self.xdata = [HC / x for x in self.xdata]
                self.xunits = SpectralUnits.KEV
            else:
                pass
        elif from_units is SpectralUnits.ANGSTROM:
            if to_units is SpectralUnits.WAVENUMBER:
                self.xdata = [1e8 / x for x in self.xdata]
                self.xunits = SpectralUnits.WAVENUMBER
            elif to_units is SpectralUnits.MICRON:
                self.xdata = [x / 1e4 for x in self.xdata]
                self.xunits = SpectralUnits.MICRON
            elif to_units is SpectralUnits.NANOMETER:
                self.xdata = [x / 1e1 for x in self.xdata]
                self.xunits = SpectralUnits.NANOMETER
            elif to_units is SpectralUnits.EV:
                self.xdata = [HC * 1e4 / x for x in self.xdata]
                self.xunits = SpectralUnits.EV
            elif to_units is SpectralUnits.KEV:
                self.xdata = [HC / x / 1e4 for x in self.xdata]
                self.xunits = SpectralUnits.KEV
            else:
                pass
        elif from_units is SpectralUnits.EV:
            if to_units is SpectralUnits.WAVENUMBER:
                self.xdata = [1e4 / x / HC for x in self.xdata]
                self.xunits = SpectralUnits.WAVENUMBER
            elif to_units is SpectralUnits.NANOMETER:
                self.xdata = [1e3 * x / HC for x in self.xdata]
                self.xunits = SpectralUnits.NANOMETER
            elif to_units is SpectralUnits.MICRON:
                self.xdata = [x / HC for x in self.xdata]
                self.xunits = SpectralUnits.MICRON
            elif to_units is SpectralUnits.ANGSTROM:
                self.xdata = [1e4 * x / HC for x in self.xdata]
                self.xunits = SpectralUnits.ANGSTROM
            elif to_units is SpectralUnits.KEV:
                self.xdata = [x / 1e3 for x in self.xdata]
                self.xunits = SpectralUnits.KEV
            else:
                pass
        elif from_units is SpectralUnits.KEV:
            if to_units is SpectralUnits.WAVENUMBER:
                self.xdata = [1e1 / x / HC for x in self.xdata]
                self.xunits = SpectralUnits.WAVENUMBER
            elif to_units is SpectralUnits.NANOMETER:
                self.xdata = [1e6 * x / HC for x in self.xdata]
                self.xunits = SpectralUnits.NANOMETER
            elif to_units is SpectralUnits.MICRON:
                self.xdata = [x * 1e3 / HC for x in self.xdata]
                self.xunits = SpectralUnits.MICRON
            elif to_units is SpectralUnits.ANGSTROM:
                self.xdata = [1e7 * x / HC for x in self.xdata]
                self.xunits = SpectralUnits.ANGSTROM
            elif to_units is SpectralUnits.EV:
                self.xdata = [x * 1e-3 for x in self.xdata]
                self.xunits = SpectralUnits.EV
            else:
                pass


if __name__ == '__main__':
    args = parse_cmdln_args()
