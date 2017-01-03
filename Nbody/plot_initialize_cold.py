#This macro plots the results in a 2d heatmap: K vs. migration rate, with the color being the various amplitudes (period ratio, eccentricity, etc).

import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
from progress.bar import Bar

#Specify directory
dir = sys.argv[1]
files = glob.glob(dir+'*.txt')
N = len(files)
i=0
while i < N:            #just want the main .txt files
    f = files[i]
    string = f.split("_")
    if string[-1]=="info.txt" or string[-1]=="elapsedtime.txt" or string[-2]=="eiasnapshot":
        files.remove(files[i])
        N -= 1
    else:
        i += 1

#figures
size=25
colorbar = 'autumn'
fig, axes = plt.subplots(nrows=3, ncols=1, figsize=(10,10), sharex=True)
plt.subplots_adjust(hspace = 0)

bar = Bar('Processing', max=N)

endpoints = 1000
A_P = []    #period amplitude
A_e1 = []   #eccentricity amplitude of planet 1
A_e2 = []   #eccentricity amplitude of planet 2
MR = []     #migration rate
K = []      #Lee & Peale (2002)
for i,f in enumerate(files):
    fos = open(f, 'r')
    time, dE, N, mig_rate, dampratio, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
    if np.any(e2>1) or np.any(e2<0)==True or np.any(e1>1) or np.any(e1<0):
        print "\nerror in %s, skipping file."%f
    else:
        P = (a2[-endpoints:]/a1[-endpoints:])**1.5
        A_P.append((np.max(P) - np.min(P))/2.)
        A_e1.append((np.max(e1[-endpoints:]) - np.min(e1[-endpoints:]))/2.)
        A_e2.append((np.max(e2[-endpoints:]) - np.min(e2[-endpoints:]))/2.)
        MR.append(mig_rate[0])
        K.append(dampratio[0])
    bar.next()

bar.finish()
sc = axes[0].scatter(K, MR, c=A_P, s=size, cmap=colorbar, lw=0)
sc2 = axes[1].scatter(K, MR, c=A_e1, s=size, cmap=colorbar, lw=0)
sc3 = axes[2].scatter(K, MR, c=A_e2, s=size, cmap=colorbar, lw=0)
axes[0].set_xscale('log')
axes[0].set_xlim([1e-2,1e6])
axes[0].set_yscale('log')
axes[1].set_yscale('log')
axes[2].set_yscale('log')
axes[0].set_ylabel('K')
axes[1].set_ylabel('K')
axes[2].set_ylabel('K')
axes[2].set_xlabel('Outer Planet Migration Speed')
plt.colorbar(sc, ax=axes[0], label='Period Ratio Amplitude')
plt.colorbar(sc2, ax=axes[1], label='e1 Amplitude')
plt.colorbar(sc3, ax=axes[2], label='e2 Amplitude')
axes[0].set_title('Variable Amplitudes: V$_{max}$ - V$_{mean}$')

plt.savefig(dir+"amplitude.png")
plt.show()

'''
    sc = axes[0].plot(MR, A_P, '.', ms=size)#, cmap=colorbar)
    #plt.colorbar(sc, ax=axes[0])
    axes[1].plot(MR, A_e1, '.', ms=size, label='e1')
    axes[1].plot(MR, A_e2, '.', ms=size, label='e2')
    axes[0].set_ylabel('Period Ratio amplitude (Prat.max - Prat.mean)')
    axes[0].set_title(r'$\tau_a/\tau_e=%.0f$'%K[0])
    axes[1].set_ylabel('e amplitude (e.max - e.mean)')
    axes[1].set_xscale('log')
    axes[1].set_xlim([1e2,1e7])
    axes[1].set_xlabel('Migration Rate')
    axes[1].legend(loc='upper left', numpoints=1)
    '''
