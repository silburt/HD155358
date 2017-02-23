# This script plots the distribution of saved_output/roundX/good_ones/*.txt, i.e. the best fits from the simulated migration runs.

import matplotlib.pyplot as plt
import numpy as np
import sys
import glob
from progress.bar import Bar

dir = sys.argv[1]
files = glob.glob(dir+'*.txt')

tau_a2, tau_e1, tau_e2 = [], [], []

bar = Bar('Processing', max=len(files))
for f in files:
    time, dE, N, mig_rate, damp1, damp2, migtime, DT, a1, e1, a2, e2, phi1, phi2, phi3, m1, m2, taua1, taue1, taua2, taue2 = np.loadtxt(open(f, 'r'), delimiter=',', unpack=True)
    tau_a2.append(taua2[0]), tau_e1.append(taue1[0]), tau_e2.append(taue2[0])
    bar.next()
bar.finish()

size=50
colorbar = 'winter'
fontsize=16
sc = plt.scatter(tau_e1,tau_e2,c=np.log10(np.asarray(tau_a2)), s=size, cmap=colorbar, lw=0, label='simulated points')
plt.colorbar(sc, label=r'log10($\tau_{a_2}$)')
plt.xlabel(r'$\tau_{e_1}$',fontsize=fontsize)
plt.ylabel(r'$\tau_{e_2}$',fontsize=fontsize)
plt.yscale('log')
plt.xscale('log')
plt.ylim([min(tau_e2)/2, max(tau_e2)*2])
plt.savefig("%splot_good_ones.pdf"%dir)
