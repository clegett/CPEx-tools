#!/usr/bin/env python3
"""continuum_removal.py: takes x and y data, a startx and stopx, and removes
the continuum of the data between startx and stopx."""

import argparse
import sys
import itertools
import os
import csv

__author__ = 'Carey Legett'
__copyright__ = 'Copyright 2018, Stony Brook University'
__credits__ = ['Carey Legett']

__version__ = '1.0'
__maintainer__ = 'Carey Legett'
__email__ = 'carey.legett@stonybrook.edu'
__status__ = 'Development'

parser = argparse.ArgumentParser(prog='continuum_removal.py',
        description='''This program takes a delimited text file containing\
        x and y values and fits a line from y(xstart) to y(xstop) and\
        normalizes all y data such that the linear fit line is at a value\
        of 1. It the writes an output file as a CSV containing "x, [y,]\
        [continuum,] cont_removed_y" data.''', epilog='''For more information\
        contact Chip at carey.legett@stonybrook.edu.''')
parser.add_argument('-f', '--input-file', required=True, help='''The\
        name of the input file (include path if the file is not in this\
        directory).''', dest='ifile')
parser.add_argument('-a', '--xstart', type=int, required=True,
        help='The x value at which to begin the continuum calculation.', 
        dest='xa')
parser.add_argument('-b', '--xstop', type=int, required=True,
        help='The x value at which to end the continuum calculation.', 
        dest='xb')
group = parser.add_mutually_exclusive_group(required=True)
group.add_argument('-d', '--delimiter', help='''The\
        delimiter for the input file. A single character representing the\ 
        actual delimiter. Mutually exclusive with -w or --whitespace''', 
        dest='delim')
group.add_argument('-w', '--whitespace', help='''Indicates the delimiter is \
        either tabs or spaces. Mutually exclusive with -d or --delimiter.''',
        action='store_true', dest='whitespace')
parser.add_argument('-n', '--skip-lines', default=0, help='''Skip\
        this many lines of the input file before reading in data.''', type=int,
        dest='skip')
parser.add_argument('-o', '--output-file', nargs='?', default='output.txt', 
        help='''The name of the output file (include path if you want the file\
        somewhere other than in the current directory). Defaults to \
        output.txt.''', dest='ofile')
parser.add_argument('-i', '--include-input', action='store_true', dest='i',
        help='Write the original input data before the (optional) continuum \
        fit and continuum removed data.')
parser.add_argument('-c', '--continuum-output', action='store_true', dest='c',
        help='Write the continuum fit for each input x value after the\
        (optional) original input data and before the continuum removed data.')
parser.add_argument('-v', '--verbose', action='store_true', dest='v', 
        help='Give status information during processing.')

args = parser.parse_args()

if args.v: print('Reading input file...')

if sys.version_info[0] == 2:
    raccess = 'rb'
    waccess = 'wb'
    kwargs = {}
else:
    raccess = 'rt'
    waccess = 'w'
    kwargs = {'newline':''}

try:
    with open(args.ifile, raccess, **kwargs) as infile:
        if args.whitespace:
            rawdata=[line.split() for line in itertools.islice(infile,args.skip,
                None)]
        else:
            rawdata=[line.split(args.delim) for line in \
                itertools.islice(infile, args.skip, None)]
    if not rawdata:
        sys.exit('No data in file ' + infile)
except IOError as e:
    sys.exit('I/O error: file {}: {}'.format(args.ifile, e))
except:
    sys.exit('Unexpected error: {}'.format(sys.exc_info()))

try:
    data = [[float(astring) for astring in inner] for inner in rawdata]
except:
    sys.exit('''Converting strings to floats failed. Check input data and\
    delimiter. Error: {}'''.format(sys.exc_info()))

if args.v: print('Calculating continuum line')

startx=float(args.xa)
stopx=float(args.xb)
starty=data[args.xa-int(data[0][0])][1]
stopy=data[args.xb-int(data[0][0])][1]

if args.v: 
    print('Start: (' + str(startx) + ',' + str(starty) + ')')
    print('Stop: (' + str(stopx) + ',' + str(stopy) + ')')

slope=(stopy-starty)/(stopx-startx)
intercept=starty-(slope*startx)

if args.v:
    print('Continuum equation: y=' + str(slope) + '*x+' + str(intercept))
    print('Removing continuum')

output = []
for element in data:
    continuum=slope*element[0]+intercept
    tempoutput = [element[0]]
    if args.i:
        tempoutput.append(element[1])
    if args.c: 
        tempoutput.append(continuum)
    tempoutput.append(element[1]/continuum)
    output.append(tempoutput)

if args.v: print('Writing output file')

try:
    with open(args.ofile, waccess, **kwargs) as outfile:
        writer = csv.writer(outfile)
        header = ['x']
        if args.i:
            header.append('original')
        if args.c:
            header.append('continuum')
        header.append('cont_removed')
        writer.writerow(header)
        writer.writerows(output)
except:
   sys.exit('Output file error: {}'.format(sys.exc_info()))
