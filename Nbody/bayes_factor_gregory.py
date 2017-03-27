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
    samples = np.load(f.split('_lnprob.npy')[0]+'_chain.npy')[:,500:,:].reshape((-1, n_params))
    lnL = np.load(f)[:,500:].reshape((-1))
    MLE = np.exp(lnL[np.argmax(lnL)])
    #std = 1
    #for i in range(n_params):
    #    std *= np.std(samples[:,i])
    std = np.linalg.det(np.cov(samples, rowvar=0)) * (2*np.pi)**(n_params/2)
    glob_prob = np.append(glob_prob, MLE*std)
    name.append(f.split('/')[-1].split('_lnprob')[0].split('taueinner_')[1])

x=range(1,len(glob_prob)+1)
ref_value = 45               #Find X largest Bayes' Factor, to be used as the reference model.
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
    #print f, taua2[0]/(2*np.pi), taue1[0]/(2*np.pi), taue2[0]/(2*np.pi)

#plotting
size=50
colorbar = 'winter'
fontsize=16
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8,10))
plt.subplots_adjust(hspace = 0.2)
sc = ax[0].scatter(K_1,K_2,c=np.log10(np.asarray(tau_a2)), s=size, cmap=colorbar, lw=0, label='simulated points')
ax[0].set_xlim([min(K_1)/2, max(K_1)*2])
ax[0].set_ylim([min(K_2)/2, max(K_2)*2])
ax[0].set_xlabel(r'$K_1$',fontsize=fontsize)
ax[0].set_ylabel(r'$K_2$',fontsize=fontsize)
ax[0].set_yscale('log')
ax[0].set_xscale('log')

sc = ax[1].scatter(tau_e1,tau_e2,c=np.log10(np.asarray(tau_a2)), s=size, cmap=colorbar, lw=0, label='simulated points')
ax[1].set_xlim([min(tau_e1)/2, max(tau_e1)*2])
ax[1].set_ylim([min(tau_e2)/2, max(tau_e2)*2])
ax[1].set_xlabel(r'$\tau_{e_1}$',fontsize=fontsize)
ax[1].set_ylabel(r'$\tau_{e_2}$',fontsize=fontsize)
ax[1].set_yscale('log')
ax[1].set_xscale('log')

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(sc, label=r'log10($\tau_{a_2}$)',cax=cbar_ax)
plt.savefig('%splot_good_ones.pdf'%dir)
