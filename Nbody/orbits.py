#This checks whether the planets are still in resonance or not. 
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
time, dE, N, mig_rate, dampratio, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
#time, dE, N, N_mini, HSF, m1, m2, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
#time, dE, N, N_mini, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)

ms=3
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,10), sharex=True)
plt.subplots_adjust(hspace = 0.35)

#plot
axes[0].plot(time, dE,'o',ms=ms, markeredgecolor='none')
axes[1].plot(time, a1, 'o', ms=ms, markeredgecolor='none')
axes[1].plot(time, a2, 'o', ms=ms, markeredgecolor='none')
#axes[1].plot(time, N_mini, 'o', ms=ms, markeredgecolor='none')
axes[2].plot(time, (a2/a1)**(1.5), 'o', ms=ms, markeredgecolor='none')
axes[3].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='planet 1')
axes[3].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='planet 2')

#labelling
axes[0].set_ylabel('Fractional Energy', fontsize=13)
axes[0].set_yscale('log')
axes[1].set_ylabel('Semi-major Axis', fontsize=13)
axes[1].set_xlabel('Time (Years)', fontsize=13)
axes[2].set_ylabel('Period Ratio', fontsize=13)
axes[2].set_xlabel('Time (Years)', fontsize=13)
axes[3].set_ylabel('Eccentricity', fontsize=13)
axes[3].set_xlabel('Time (Years)', fontsize=13)

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'_orbit.png')
#plt.show()
