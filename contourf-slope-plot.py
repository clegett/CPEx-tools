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
z = np.array([3.20956284279089E-05, 5.45307953306509E-05, 3.163518322858E-05,
        4.81510400808949E-05, 6.10446151362146E-05, 3.09635744172054E-05,
        4.83242275634782E-05, 4.35584011140235E-05, 5.42038128706957E-05,
        4.95945627071086E-05, 3.28552018772379E-06, 5.24805527199726E-05,
        4.17476132997527E-05, 2.18048627359877E-06])

xi = np.linspace(min(x), max(x), 100)
yi = np.linspace(min(y), max(y), 100)
zi = griddata(x, y, z, xi, yi, interp='linear')

fig = plt.figure(facecolor='black')
font = {'size' : 20}
plt.rc('font', **font)
CS = plt.contourf(xi, yi, zi, 60, vmax=0.00007, vmin=0,
        cmap=plt.cm.plasma)
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
cbar.ax.set_ylabel('800-1700 nm linear fit slope', color='white')
cbar.outline.set_edgecolor('white')
cbar.ax.yaxis.set_tick_params(color='black', labelsize=0)
plt.setp(plt.getp(cbar.ax.axes, 'yticklabels'), color='black')
plt.show()
