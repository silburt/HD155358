#This checks whether the planets are still in resonance or not. 
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])

fos = open(''+file_name, 'r')
time, dE, a1, e1, l1, w1, M1, a2, e2, l2, w2, M2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)

ms=6
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,13), sharex=True)
plt.subplots_adjust(hspace = 0.2)

#period
P1 = (a1[0]**3 /0.92)**0.5 * 2*np.pi
time /= P1

#plot
axes[0].plot(time, np.abs(dE),'o',ms=ms, markeredgecolor='none')
axes[1].plot(time, phi1, '.', ms=ms, markeredgecolor='none', label='phi 1')
axes[1].plot(time, phi2, '.', ms=ms, markeredgecolor='none', label='phi 2')
#axes[2].plot(time, (a2/a1)**(1.5), 'o', ms=ms, markeredgecolor='none')
axes[2].plot(time, a1, 'o', ms=ms, markeredgecolor='none', label='planet 1')
axes[2].plot(time, a2, 'o', ms=ms, markeredgecolor='none', label='planet 2')
axes[3].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='planet 1')
axes[3].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='planet 2')

#labelling
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].set_ylabel('dE', fontsize=13)
axes[1].set_ylabel('Resonant Angle', fontsize=13)
axes[1].set_xscale('log')
axes[1].legend(loc='upper left')
axes[2].set_ylabel('Semi-Major Axis', fontsize=13)
axes[2].set_xscale('log')
axes[3].set_ylabel('Eccentricity', fontsize=13)
axes[3].set_xscale('log')

name = file_name.split('.csv')[0]
plt.savefig('%s_orbit.png'%name)
#plt.show()
