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
#time, dE, N, N_mini, HSF, m1, m2, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
time, dE, N, N_mini, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)

ms=3
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,10))
plt.subplots_adjust(hspace = 0.2)

#periods
P1 = np.sqrt(4*np.pi**2 * a1**3)
P2 = np.sqrt(4*np.pi**2 * a2**3)

#plot
axes[2].plot(time,P2/P1,'o',ms=ms, markeredgecolor='none')
axes[1].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='inner')
axes[1].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='outer')
axes[0].plot(time,phi1, 'o', ms=ms, markeredgecolor='none', label='phi1')
axes[0].plot(time,phi2, 'o', ms=ms, markeredgecolor='none', label='phi2')
axes[0].plot(time,phi3, 'o', ms=ms, markeredgecolor='none', label='phi3')

#im = axes[2].scatter(e1*np.cos(phi1),e1*np.sin(phi1),c=time, cmap=cm.rainbow,lw=0)
#axes[3].scatter(e1*np.cos(phi2),e1*np.sin(phi2),c=time, cmap=cm.rainbow,lw=0)
#plt.colorbar(im, label='elapsed time (yr)')

#labelling
axes[0].set_ylabel('Resonant Angles', fontsize=13)
axes[0].legend(loc='upper left')
axes[1].set_ylabel('Eccentricity', fontsize=13)
axes[1].legend(loc='upper left')
axes[2].set_ylabel('Period Ratio', fontsize=13)
axes[2].set_xlabel('Time (Years)', fontsize=13)
#axes[2].set_ylabel('e1*cos(phi1)', fontsize=13)
#axes[2].set_xlabel('e1*sin(phi1)', fontsize=13)
#axes[3].set_ylabel('e1*cos(phi2)', fontsize=13)
#axes[3].set_xlabel('e1*sin(phi2)', fontsize=13)

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'_rescheck.png')
#plt.show()
