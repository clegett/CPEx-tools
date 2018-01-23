#!/usr/bin/env python3
"""continuum_removal.py: takes x and y data, a startx and stopx, and removes
the continuum of the data between startx and stopx."""

import argparse
import csv
import sys
import itertools

__author__ = "Carey Legett"
__copyright__ = "Copyright 2018, Stony Brook University"
__credits__ = ["Carey Legett"]

__version__ = "1.0"
__maintainer__ = "Carey Legett"
__email__ = "carey.legett@stonybrook.edu"
__status__ = "Development"

parser = argparse.ArgumentParser(prog="continuum_removal.py",
            description='''This program takes a delimited text file containing
            x and y values and fits a line from y(xstart) to y(xstop) and
            normalizes all y data such that the linear fit line is at a value
            of 1.''', epilog='''For more information, contact Chip at
            carey.legett@stonybrook.edu.''')
parser.add_argument('-f', '--input-file', nargs=1, required=True, help='''The 
    name of the input file (include path if the file is not in this directory
    )''', dest='ifile')
parser.add_argument('-o', '--output-file', nargs='?', default='output.txt', 
    help='''The name of the output file (include path if you want the file 
    somewhere other than in the current directory). Defaults to output.txt''',
    dest='ofile')
parser.add_argument('-a', '--xstart', type=float, nargs=1, required=True,
    help='The x value at which to begin the continuum calculation.', dest='xa')
parser.add_argument('-b', '--xstop', nargs=1, type=float, required=True,
    help='The x value at which to end the continuum calculation.', dest='xb')
parser.add_argument('-d', '--delimiter', nargs='?', default=',', help='''The
delimiter for the input and output files. Defaults to \',\' ''', dest='delim')
parser.add_argument('-n', '--start-line', nargs='?', default=0, help='''Skip
    this many lines of the input file before reading in data.''', type=int,
    dest='skip')

args = parser.parse_args()
print(args.delim)
mydialect = csv.register_dialect('asd', delimiter=args.delim, skipinitialspace=False)
with open(args.ifile[0], newline='')  as f:
    reader = csv.reader(itertools.islice(f,args.skip,None), dialect='asd')
    try:
        data=list(reader)
    except csv.Error as e:
        sys.exit('file {}, line {}: {}'.format(filename, reader.line_num, e))

