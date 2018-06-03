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
    start_wl = 700 * conv.ureg['nanometers']
    stop_wl = 1600 * conv.ureg['nanometers']
    start_wn = start_wl.to(1 / conv.ureg['centimeters'], 'sp').magnitude
    stop_wn = stop_wl.to(1 / conv.ureg['centimeters'], 'sp').magnitude

    temp_data = cr.remove_continuum(temp_data, start_wn, stop_wn, False)
    wls = [row[0] for row in temp_data]
    cr_spec = [row[3] for row in temp_data]
    file_data.append(list(zip(wls, cr_spec)))

# calculate band depth
band_depths = []
band_index = 400
# print(f'length of file_data: {len(file_data)}')
for run in file_data:
    # print(f'run length: {len(run)}')
    # print(f'length of run[0]: {len(run[0])}')
    # print(f'run[0]: {run[0]}')
    band_depths.append(run[band_index][1])

# write band depths to file
try:
    with open('band_depths.csv', 'w', newline='') as outfile:
        writer = csv.writer(outfile)
        header = ['file', f'banddepth_at_{band_index}']

        writer.writerow(header)
        writer.writerows(zip(filelist, [1 - depth for depth in band_depths]))
except IOError as e:
    sys.exit(f'Output file error: {e}, {sys.exc_info()}')
