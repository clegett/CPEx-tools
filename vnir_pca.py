#!/usr/bin/env python3

from sklearn.decomposition import PCA
from sklearn.preprocessing import StandardScaler
from scipy.signal import savgol_filter as sf
import pandas as pd
import csv
import sys
import matplotlib
matplotlib.use('TkAgg')
import matplotlib.pyplot as plt

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
        sys.exit('I/O error: file {}: {}'.format(filename, e))


pca = PCA(n_components=4)
infile = 'temple_matlab_format.csv'
spectral_data = read_csv_to_list(infile)

run_names = [run_map[item] for item in spectral_data[0][1:]]
wls = [row[0] for row in spectral_data[1:]]
refl = [row[1:] for row in spectral_data[1:]]
# translate the data so that each row is a sample
reflr = list(map(list, zip(*refl)))
filt_ref = [sf(row, 5, 3) for row in reflr]
filt_dref = [sf(row, 5, 3, 1) for row in reflr]
df = pd.DataFrame({'x': wls, 'reflr': reflr[0], 'filt_ref': filt_ref[0],
                   'filt_dref': filt_dref[0]})

plt.plot('x', 'reflr', data=df)
plt.plot('x', 'filt_ref', data=df)
plt.plot('x', 'filt_dref', data=df)
plt.legend()
print('trying to plot')
plt.show()

data = pd.DataFrame(dict(zip(wls, refl)), index=run_names)
data.drop(['Bad1', 'Bad2', 'Bad3'])

x = data.loc[:, wls].values
x = StandardScaler().fit_transform(x)
pc = pca.fit_transform(x)
principalDf = pd.DataFrame(data=pc)
