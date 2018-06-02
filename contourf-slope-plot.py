#!/usr/bin/env python3

import numpy as np
from matplotlib import pyplot as plt
from matplotlib.mlab import griddata

title = '25x cluster 800-1700 slope'

# x dFe
# y wt%Fe
# z linear slope

x = np.array([10, 10, 20, 20, 20, 40, 40, 40, 64, 64, 64, 86, 86, 86])
y = np.array([0.010233975665857, 0.102271913273833, 0.010233975665857,
        0.102271913273833, 1.015983940899, 0.010233975665857,
        0.102271913273833, 1.015983940899, 0.08379231030076, 0.833396601843304,
        7.9068349660572, 0.101642123980333, 1.00976847323586,
        9.47736424395359])
z = np.array([0.000145383, 0.000254705, 0.000143673, 0.000266494, 
	0.000472658, 0.000139502, 0.000269381, 0.000495233, 
	0.000300413, 0.000531951, 0.000150082, 0.000293789, 
	0.000518836, 0.00010917])

xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
zi = griddata(x, y, z, xi, yi, interp='linear')

fig = plt.figure(facecolor='black', figsize=(16, 9))
font = {'size' : 20}
plt.rc('font', **font)
CS = plt.contourf(xi, yi, zi, 60, vmax=0.00053, vmin=0.0001,
        cmap=plt.cm.plasma)
CS2 = plt.contour(CS, levels=CS.levels[::4], colors='black')
CS.ax.set_yscale('log')
CS.ax.set_facecolor((0, 0, 0))
CS.ax.spines['bottom'].set_color('white')
CS.ax.spines['top'].set_color('white')
CS.ax.spines['left'].set_color('white')
CS.ax.spines['right'].set_color('white')
for t in CS.ax.xaxis.get_ticklines(): t.set_color('white')
for t in CS.ax.yaxis.get_ticklines(): t.set_color('white')
CS.ax.xaxis.label.set_color('white')
CS.ax.yaxis.label.set_color('white')
CS.ax.tick_params(colors='white')
#CS.ax.grid(c=(0.9, 0.9, 0.9), lw=0.4, linestyle='dotted', axis='y')
#CS.ax.grid(c=(0.9, 0.9, 0.9), lw=0.2, linestyle='dotted', axis='y',
#    which='minor')
plt.scatter(x, y, marker='s', s=2**7, c='w', zorder=10, edgecolors='black')
plt.title('{}'.format(title))
plt.xlabel('Fe particle diameter (nm)')
plt.ylabel('wt% Fe')
cbar = plt.colorbar(CS)
cbar.ax.set_ylabel('linear fit slope', color='white')
cbar.outline.set_edgecolor('white')
cbar.ax.yaxis.set_tick_params(color='black', labelsize=0)
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='black')
plt.show()
fig.savefig('contourslope.png', edgecolor='black', facecolor='black')
