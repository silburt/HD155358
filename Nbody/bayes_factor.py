#Once you've isolated the good runs, it's time to compute the Bayes' factor for each run, i.e. sum the likelihood values for each MCMC chain. 
#Then plot_good_ones.py has been rolled into this script now, and I plot the runs with the best Bayes' factor.

#Maybe need to generalize this routine a bit more so that you can multiply by 1/prior after, and automatically be able to load the 'prior constant', i.e. 1/(e_max - e_min) * 1/(a_max - a_min) * ...

import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
import rebound

def get_X_largest(array, X):
    min = array[np.argmin(array)]
    for i in range(X-1):
        array[np.argmax(array)] = min
    return np.argmax(array)

#initial args/settings
dir = sys.argv[1]
files = glob.glob(dir+'*_lnprob.npy')
n_params = 6  #x_s, x_t, y_s, y_t, phi, jitter2
glob_prob = np.empty(0)
name = []

#main loop
print "Getting Bayes' Factors..."
for f in files:
    lnL = np.load(f)[:,500:].reshape((-1))
    glob_prob = np.append(glob_prob, np.sum(np.exp(lnL)))
    name.append(f.split('/')[-1].split('_lnprob')[0].split('taueinner_')[1])

x=range(1,len(glob_prob)+1)
ref_value = 6               #Find X largest Bayes' Factor, to be used as the reference model.
bayes_factor = glob_prob/glob_prob[get_X_largest(np.copy(glob_prob),ref_value)]

plt.plot(x, bayes_factor, 'o')
plt.xticks(x, name, rotation='vertical')
plt.xlim([0,x[-1]+1])
plt.yscale('log')
plt.ylabel('posterior odds ratio (normalized to %d highest)'%ref_value)
plt.savefig('%sbayes_factor.pdf'%dir)
plt.clf()

###################################
#plot values with bayes_factor > 1
print "Plotting runs with good Bayes' Factors..."
good_files = np.asarray(files)[bayes_factor >=1]

tau_a2, tau_e1, tau_e2, K_1, K_2 = [], [], [], [], []
for f in good_files:
    time, dE, N, mig_rate, damp1, damp2, migtime, DT, a1, e1, a2, e2, phi1, phi2, phi3, m1, m2, taua1, taue1, taua2, taue2 = np.loadtxt(open(f.split('_lnprob.npy')[0]+'.txt', 'r'), delimiter=',', unpack=True)
    tau_a2.append(taua2[0]/(2*np.pi)), tau_e1.append(taue1[0]/(2*np.pi)), tau_e2.append(taue2[0]/(2*np.pi)), K_1.append(damp1[0]), K_2.append(damp2[0])
    print f, taua2[0]/(2*np.pi), damp1[0], damp2[0]

size=50
colorbar = 'winter'
fontsize=16
sc = plt.scatter(K_1,K_2,c=np.log10(np.asarray(tau_a2)), s=size, cmap=colorbar, lw=0, label='simulated points')
plt.colorbar(sc, label=r'log10($\tau_{a_2}$)')
plt.xlabel(r'$K_1$',fontsize=fontsize)
plt.ylabel(r'$K_2$',fontsize=fontsize)
plt.yscale('log')
plt.xscale('log')
plt.ylim([min(K_2)/2, max(K_2)*2])
plt.savefig('%splot_good_ones.pdf'%dir)

#sc = plt.scatter(tau_e1,tau_e2,c=np.log10(np.asarray(tau_a2)), s=size, cmap=colorbar, lw=0, label='simulated points')
#plt.ylim([min(tau_e2)/2, max(tau_e2)*2])
