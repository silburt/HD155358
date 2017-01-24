#This checks whether the planets are still in resonance or not. 
import sys
import matplotlib.pyplot as plt
import numpy as np
import math
import matplotlib.cm as cm
import re

colors=['b','g','m','r','c','y']

file_name=str(sys.argv[1])
time, dE, RV = np.loadtxt(open(file_name, 'r'), delimiter=',', unpack=True)

ms=3
fig, axes = plt.subplots(nrows=2, ncols=1, figsize=(10,10), sharex=True)
plt.subplots_adjust(hspace = 0.35)

#plot
axes[0].plot(time, dE,'o',ms=ms, markeredgecolor='none')
axes[1].plot(time, RV, ms=ms, markeredgecolor='none')

#labelling
axes[0].set_ylabel('Fractional Energy', fontsize=13)
axes[0].set_yscale('log')
axes[1].set_ylabel('Radial Velocity', fontsize=13)
axes[1].set_xlabel('Time (days)', fontsize=13)

file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'.png')
plt.show()
