#This is an "orbits.py" plotter for the stability.py runs.
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
time, dE, a1, e1, l1, w1, M1, a2, e2, l2, w2, M2 = np.loadtxt(fos, delimiter=',', unpack=True)

ms=6
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(8,10), sharex=True)
plt.subplots_adjust(hspace = 0.2)

P1 = (a1[0]**3 /0.92)**0.5 * 2*np.pi
time /= P1

#plot
axes[0].plot(time, np.abs(dE), '.', ms=ms, markeredgecolor='none', label='phi 1')
axes[1].plot(time, a1, 'o', ms=ms, markeredgecolor='none')
axes[1].plot(time, a2, 'o', ms=ms, markeredgecolor='none')
axes[2].plot(time, (a2/a1)**(1.5), 'o', ms=ms, markeredgecolor='none')
axes[3].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='planet 1')
axes[3].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='planet 2')

axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[1].set_xscale('log')
axes[2].set_xscale('log')
axes[3].set_xscale('log')

axes[0].set_ylabel('Fractional Energy Error')
axes[1].set_ylabel('a')
axes[2].set_ylabel('Period Ratio')
axes[3].set_ylabel('e')
axes[3].set_xlabel('time (years)')

plt.show()
