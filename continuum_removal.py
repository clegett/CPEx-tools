#!/usr/bin/env python3
"""continuum_removal.py: takes x and y data, a startx and stopx, and removes
the continuum of the data between startx and stopx."""

import argparse
import sys
import itertools
import os

__author__ = 'Carey Legett'
__copyright__ = 'Copyright 2018, Stony Brook University'
__credits__ = ['Carey Legett']

__version__ = '1.0'
__maintainer__ = 'Carey Legett'
__email__ = 'carey.legett@stonybrook.edu'
__status__ = 'Development'

parser = argparse.ArgumentParser(prog='continuum_removal.py',
        description='''This program takes a delimited text file containing
        x and y values and fits a line from y(xstart) to y(xstop) and
        normalizes all y data such that the linear fit line is at a value
        of 1.''', epilog='''For more information, contact Chip at
        carey.legett@stonybrook.edu.''')
parser.add_argument('-f', '--input-file', required=True, help='''The
        name of the input file (include path if the file is not in this
        directory)''', dest='ifile')
parser.add_argument('-o', '--output-file', nargs='?', default='output.txt', 
        help='''The name of the output file (include path if you want the file
        somewhere other than in the current directory). Defaults to 
        output.txt''', dest='ofile')
parser.add_argument('-a', '--xstart', type=int, required=True,
        help='The x value at which to begin the continuum calculation.', 
        dest='xa')
parser.add_argument('-b', '--xstop', type=int, required=True,
        help='The x value at which to end the continuum calculation.', 
        dest='xb')
parser.add_argument('-d', '--delimiter', nargs='?', default=',', help='''The
        delimiter for the input and output files. Defaults to \',\' ''', 
        dest='delim')
parser.add_argument('-n', '--start-line', nargs='?', default=0, help='''Skip
        this many lines of the input file before reading in data.''', type=int,
        dest='skip')

args = parser.parse_args()

try:
    with open(args.ifile, newline='')  as f:
        rawdata=itertools.islice(f.readlines(),args.skip,None)
        f.close()
    if not rawdata:
        sys.exit('No data in file ' + f)
except IOError as e:
    sys.exit('I/O error: file {}: {}'.format(args.ifile[0], e))
except:
    sys.exit('Unexpected error: {}'.format(sys.exc_info()[0]))

strippeddata = map(str.strip, rawdata)
strdata = [line.split() for line in strippeddata]

try:
    data = [[float(astring) for astring in inner] for inner in strdata]
except:
    sys.exit('''Converting strings to floats failed. Check input data. Error:
    {}'''.format(sys.exc_info()[0]))


startx=float(args.xa)
stopx=float(args.xb)
starty=data[args.xa-int(data[0][0])][1]
stopy=data[args.xb-int(data[0][0])][1]

slope=(stopy-starty)/(stopx-startx)
intercept=starty-(slope*startx)

newdata = list()

for i in range(len(data)):
    continuum=slope*data[i][1]+intercept
    newdata[i][0]=data[i][0]
    newdata[i][1]=data[i][1]/continuum
    print(data[i] + newdata[i])
