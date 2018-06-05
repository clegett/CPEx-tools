#!/usr/bin/env python3

import sys
import csv
import plotly.plotly as py
import plotly.graph_objs as go
from scipy.stats import linregress as lin

run_map = {'Avg-1': 'Ol', 'Avg-2': 'Ol05', 'Avg-3': 'Ol1', 'Avg-4': 'Ol2',
           'Avg-5': 'Ol5', 'Avg-6': 'Bad1', 'Avg-6-2': 'Px', 'Avg-7': 'Px05',
           'Avg-8': 'Px1', 'Avg-9': 'Px2', 'Avg-10': 'Bad2', 'Avg-10-2':
           'Px5', 'Avg-11': 'Bad3', 'Avg-11-2': 'Fs', 'Avg-12': 'Fs05',
           'Avg-13': 'Fs1', 'Avg-14': 'Fs2', 'Avg-15': 'Fs5', 'Avg-16':
           'Fs01Ol', 'Avg-17': 'Fs01Ol05', 'Avg-18': 'Fs01Ol1', 'Avg-19':
           'Fs01Ol2', 'Avg-20': 'Fs01Ol5', 'Avg-21': 'Fs05Ol', 'Avg-22':
           'Fs05Ol05', 'Avg-23': 'Fs05Ol1', 'Avg-24': 'Fs05Ol2', 'Avg-25':
           'Fs05Ol5', 'Avg-26': 'Fs1Ol', 'Avg-27': 'Fs1Ol05', 'Avg-28':
           'Fs1Ol1', 'Avg-29': 'Fs1Ol2', 'Avg-30': 'Fs1Ol5', 'Avg-31':
           'Fs2Ol', 'Avg-32': 'Fs2Ol05', 'Avg-33': 'Fs2Ol1', 'Avg-34':
           'Fs2Ol2', 'Avg-35': 'Fs2Ol5', 'Avg-36': 'Fs5Ol', 'Avg-37':
           'Fs5Ol05', 'Avg-38': 'Fs5Ol1', 'Avg-39': 'Fs5Ol2', 'Avg-40':
           'Fs5Ol5', 'Avg-41': 'Fs01Px', 'Avg-42': 'Fs01Px05', 'Avg-43':
           'Fs01Px1', 'Avg-44': 'Fs01Px2', 'Avg-45': 'Fs01Px5', 'Avg-46':
           'Fs05Px', 'Avg-47': 'Fs05Px05', 'Avg-48': 'Fs05Px1', 'Avg-49':
           'Fs05Px2', 'Avg-50': 'Fs05Px5', 'Avg-51': 'Fs1Px', 'Avg-52':
           'Fs1Px05', 'Avg-53': 'Fs1Px1', 'Avg-54': 'Fs1Px2', 'Avg-55':
           'Fs1Px5', 'Avg-56': 'Fs2Px', 'Avg-57': 'Fs2Px05', 'Avg-58':
           'Fs2Px1', 'Avg-59': 'Fs2Px2', 'Avg-60': 'Fs2Px5', 'Avg-61': 'Fs5Px',
           'Avg-62': 'Fs5Px05', 'Avg-63': 'Fs5Px1', 'Avg-64': 'Fs5Px2',
           'Avg-65': 'Fs5Px5'}


def read_csv_to_list(filename):
    try:
        with open(filename, 'rt', newline='') as f:
            reader = csv.reader(f)
            return list(reader)
    except IOError as e:
        sys.exit(f'I/O error: file {filename}: {e}')


def write_list_to_csv(filename, data, header=None):
    try:
        with open(filename, 'w', newline='') as o:
            writer = csv.writer(o)
            if header is not None:
                writer.writerow(header)
            writer.writerows(data)
    except IOError as e:
        sys.exit(f'I/O error: file {filename}: {e}')


def transpose_list(alist):
    return list(map(list, zip(*alist)))


def make_axis(title, tickangle):
    return {
      'title': title,
      'titlefont': {'size': 20},
      'tickangle': tickangle,
      'tickfont': {'size': 15},
      'tickcolor': 'rgba(0,0,0,0)',
      'ticklen': 5,
      'showline': True,
      'showgrid': True
    }


infile = 'temple_matlab_format.csv'
spectral_data = read_csv_to_list(infile)
run_names = [run_map[item] for item in spectral_data[0][1:]]
wls = [float(row[0]) for row in spectral_data[1:]]
wls = wls[100:-675]
# print(wls)
refl = [[float(item)for item in row[1:]] for row in spectral_data[1:]]
# translate the data so that each row is a sample
reflr = transpose_list(refl)
reflr = [sample[100:-675] for sample in reflr]
# print([len(reflr[0]), len(wls)])

stats = []
flattened = []
flat_stacked = []
for count, run in enumerate(reflr):
    this_slope = lin(wls, run)[0]
    flat_run = []
    sum_flat = 0
    for index, item in enumerate(run):
        ref = item - (this_slope * wls[index])
        sum_flat = sum_flat + ref
        flat_run.append(ref)
    avg = sum_flat / len(run)
    run_max = max(run)
    run_min = min(run)
    run_delta = run_max - run_min
    stats.append([run_names[count], this_slope, avg, run_max, run_min,
                  run_delta])
    #print(f'len run: {len(run)}, avg: {avg}')
    flat_stacked.append([wl - avg for wl in flat_run])
    flattened.append(flat_run)
data_dict = [{'slope': row[1], 'albedo': row[2], 'contrast': row[5], 'label':
              row[0]} for row in stats]
slope_values = [x['slope'] for x in data_dict]
albedo_values = [x['albedo'] for x in data_dict]
contrast_values = [x['contrast'] for x in data_dict]
max_slope = max(slope_values)
min_slope = min(slope_values)
max_albedo = max(albedo_values)
min_albedo = min(albedo_values)
max_contrast = max(contrast_values)
min_contrast = min(contrast_values)

normed_dict = [{'slope': (row['slope'] - min_slope)/(max_slope - min_slope),
                'albedo': (row['albedo'] - min_albedo)/(max_albedo -
                                                        min_albedo),
                'contrast': (row['contrast'] - min_contrast)/(max_contrast -
                                                              min_contrast),
                'label': row['label']} for row in data_dict if 'Bad' not in
               data_dict['label']]


fdata = [{
    'type': 'scatterternary',
    'mode': 'markers',
    'a': [i for i in map(lambda x: x['slope'], normed_dict)],
    'b': [i for i in map(lambda x: x['albedo'], normed_dict)],
    'c': [i for i in map(lambda x: x['contrast'], normed_dict)],
    'text': [i for i in map(lambda x: x['label'], normed_dict)],
    'marker': {
        'symbol': 100,
        'color': '#DB7365',
        'size': 14,
        'line': {'width': 2}}, }]

layout = {
    'ternary': {
        'sum': 1,
        'aaxis': make_axis('slope', 0),
        'baxis': make_axis('<br>albedo', 45),
        'caxis': make_axis('<br>contrast', -45)
    },
    'annotations': [{'showarrow': False, 'text':
                     'Simple Ternary Plot with Markers', 'x': 0.5, 'y': 1.3,
                     'font': {'size': 15}}]}

fig = {'data': fdata, 'layout': layout}
py.iplot(fig, validate=False)


#print(flattened[0])
# print(flat_stacked[0])
header = ['wavelength'] + run_names
t_flat = transpose_list(flattened)
t_flat_stack = transpose_list(flat_stacked)
for index, wl in enumerate(wls):
    t_flat[index] = [wl] + t_flat[index]
    t_flat_stack[index] = [wl] + t_flat_stack[index]
# print(t_flat[0])
write_list_to_csv('stats.csv', transpose_list(stats))
write_list_to_csv('flattened.csv', t_flat, header=header)
write_list_to_csv('flat_stacked.csv', t_flat_stack, header=header)
