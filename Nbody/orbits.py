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
time, dE, N, mig_rate, damp1, damp2, migtime, DT, a1, e1, a2, e2, phi1, phi2, phi3, m1, m2, taua1, taue1, taua2, taue2 = np.loadtxt(fos, delimiter=',', unpack=True)

ms=6
fig, axes = plt.subplots(nrows=5, ncols=1, figsize=(10,13), sharex=True)
plt.subplots_adjust(hspace = 0.2)

P1 = (a1[0]**3 /0.92)**0.5 * 2*np.pi
time /= P1

#plot
axes[0].plot(time, phi1, '.', ms=ms, markeredgecolor='none', label='phi 1')
axes[0].plot(time, phi2, '.', ms=ms, markeredgecolor='none', label='phi 2')
#axes[0].plot(time, phi3, '.', ms=ms, markeredgecolor='none', label='phi 3')
#axes[0].plot(time, dE,'o',ms=ms, markeredgecolor='none')
axes[1].plot(time, a1, 'o', ms=ms, markeredgecolor='none')
axes[1].plot(time, a2, 'o', ms=ms, markeredgecolor='none')
axes[2].plot(time, (a2/a1)**(1.5), 'o', ms=ms, markeredgecolor='none')
axes[3].plot(time, e1, 'o', ms=ms, markeredgecolor='none', label='planet 1')
axes[3].plot(time, e2, 'o', ms=ms, markeredgecolor='none', label='planet 2')

P2 = a2**1.5
axes[4].plot(time, taue1/P2, 'o', ms=2*ms, color='black', markeredgecolor='none', label='taue planet 1')
axes[4].plot(time, taue2/P2, 'o', ms=ms, color='green', markeredgecolor='none', label='taue planet 2')
axes[4].plot(time, taua2/P2, 'o', ms=ms/3, color='red', markeredgecolor='none', label='taua planet 1')
axes[4].set_yscale('log')

plot_bounds = 0
if plot_bounds == 1:
    mig_time = migtime[0]/P1
    D_T = DT[0]/P1
    axes[0].plot([mig_time, mig_time], [0, 2*np.pi], 'r--',lw=4)
    axes[0].plot([D_T, D_T], [0, 2*np.pi], 'c--',lw=4)
    axes[1].plot([mig_time, mig_time], [0, max(a2)], 'r--',lw=4)
    axes[1].plot([D_T, D_T], [0, max(a2)], 'c--',lw=4)
    axes[2].plot([mig_time, mig_time], [1, 3], 'r--',lw=4)
    axes[2].plot([D_T, D_T], [1, 3], 'c--',lw=4)
    axes[3].plot([mig_time, mig_time], [0, max([max(e1),max(e2)])], 'r--',lw=4)
    axes[3].plot([D_T, D_T], [0, max([max(e1),max(e2)])], 'c--',lw=4)

#Goldreich & Schlichting Eccentricity Theory
#a, b, c = 0.630, 1.19, 0.428
#taun1 = -2*taua1/3
#e2_th = ( taue2/(6*taun1) )**0.5 * ( (1 - m2/m1)/( (1 + m2/(2*a*m1))*(1 + (c/b)**2/(4*a)) ) )**0.5
#e1_th = e2_th/(2*(a*b/c)*m1/m2)
#axes[3].plot(time, e1_th, color='black', label='e1_theory')
#axes[3].plot(time, e2_th, 'k--', label='e2_theory')
#axes[3].legend(loc='upper left', fontsize=10)
#print np.mean(e1_th), np.mean(e2_th)

#labelling
#axes[0].set_xscale('log')
axes[0].legend(loc='upper left')
#axes[0].set_xlim([100,max(time)])
#axes[0].set_xlim([18000,19000])
axes[0].set_ylabel('Resonant Angle', fontsize=13)
axes[0].set_xlabel('Orbital Period of Inner Planet', fontsize=13)
axes[1].set_ylabel('Semi-major axis', fontsize=13)
axes[1].set_xlabel('Orbital Period of Inner Planet', fontsize=13)
axes[2].set_ylabel('Period Ratio', fontsize=13)
axes[2].set_xlabel('Orbital Period of Inner Planet', fontsize=13)
axes[3].set_ylabel('Eccentricity', fontsize=13)
axes[3].set_xlabel('Orbital Period of Inner Planet', fontsize=13)
axes[4].set_ylabel('Damping', fontsize=13)
axes[4].set_yscale('log')


file_output_name = re.sub('\.txt$', '', file_name)
plt.savefig(file_output_name+'_orbit.png')
#plt.show()
