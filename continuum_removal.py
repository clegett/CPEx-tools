#!/usr/bin/env python3
"""Removes the continuum (background) of data over a specified range.

Removes the continuum of specified x,y data from a delimited text file by 
calculating a linear fit to two points provided by the user and normalizing 
the data such that the linear fit has a value of 1. Output writes to a CSV file
with two to four columns: "x, [y,] [continuum,] continuum_removed_y" as
controlled by command line arguments.

Args
"""
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

def parse_cmdln_args():

    parser = argparse.ArgumentParser(prog='continuum_removal.py',
            description=__doc__, epilog='''For more information\
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
    
    return(parser.parse_args())


def read_input_file(input_filename, input_delimiter, whitespace_flag, 
                    lines_to_skip, verbose_flag):
    if verbose_flag: print('Reading input file...')
    
    if sys.version_info[0] == 2:
        raccess = 'rb'
        kwargs = {}
    else:
        raccess = 'rt'
        kwargs = {'newline':''}
    
    try:
        with open(input_filename, raccess, **kwargs) as infile:
            if whitespace_flag:
                rawdata=[line.split() for line in itertools.islice(infile,
                lines_to_skip, None)]
            else:
                rawdata=[line.split(input_delimiter) for line in \
                    itertools.islice(infile, lines_to_skip, None)]
        if not rawdata:
            sys.exit('No data in file ' + input_filename)
    except IOError as e:
        sys.exit('I/O error: file {}: {}'.format(input_filename, e))
    except:
        sys.exit('Unexpected error: {}'.format(sys.exc_info()))
    
    try:
        data = [[float(astring) for astring in inner] for inner in rawdata]
    except:
        sys.exit('''Converting strings to floats failed. Check input data and\
        delimiter. Error: {}'''.format(sys.exc_info()))

    return(data)

def remove_continuum(data, xa, xb, verbose_flag):
    if args.v: print('Calculating continuum line')
    
    startx=xa
    stopx=xb
    starty=data[int(xa)-int(data[0][0])][1]
    stopy=data[int(xb)-int(data[0][0])][1]
    
    if verbose_flag: 
        print('Start: (' + str(startx) + ',' + str(starty) + ')')
        print('Stop: (' + str(stopx) + ',' + str(stopy) + ')')
    
    slope=(stopy-starty)/(stopx-startx)
    intercept=starty-(slope*startx)
    
    if verbose_flag:
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

    return(output)

def write_cont_removed_output(output, output_filename, verbose_flag):
    if verbose_flag: print('Writing output file')
    
    if sys.version_info[0] == 2:
        waccess = 'wb'
        kwargs = {}
    else:
        waccess = 'w'
        kwargs = {'newline':''}

    try:
        with open(output_filename, waccess, **kwargs) as outfile:
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

if __name__ == '__main__':
    args = parse_cmdln_args()
    data = read_input_file(args.ifile, args.delim, args.whitespace, args.skip,
                           args.v)
    output = remove_continuum(data, float(args.xa), float(args.xb), args.v)
    write_cont_removed_output(output, args.ofile, args.v)
