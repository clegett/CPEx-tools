#!/usr/bin/env python3
"""Calculate band depth of a spectrum given 3 points."""

import csv
import sys
import glob
import continuum_removal as cr
import wl_to_wn as conv

# get list of asc files in directory matching pattern
filelist = glob.glob('*.asc')

file_data = []
# read file and convert to wn
for file in filelist:
    temp_data = conv.read_input_file(file, '', True, 0, 'nanometers', False)
    temp_data = [[line[0].to(1 / conv.ureg.centimeter, 'sp').magnitude,
                  line[1]] for line in temp_data]

    # continuum remove on specified points
    temp_data = cr.remove_continuum(temp_data, 700, 1600, False)
    file_data.append([temp_data[0], temp_data[3]])

# calculate band depth
band_depths = []
band_index = 550
for run in file_data:
    band_depths.append(run[band_index][1])

# write band depths to file
try:
    with open('band_depths.csv', 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        header = ['file', f'banddepth_at_{band_index}']

        writer.writerow(header)
        writer.writerows(zip(filelist, band_depths))
except IOError as e:
    sys.exit(f'Output file error: {e}, {sys.exc_info()}')
