#This macro is going to explore whether , for a given K and migration rate, the libration amplitude will significantly change for different mass draws from my MCMC posterior distribution. Hopefully not...

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

N_runs = 100   #N_runs per (K,mig_rate) combo
K = np.concatenate((np.ones(2*N_runs),10*np.ones(2*N_runs),100*np.ones(2*N_runs)))
mig_rate = np.concatenate((320*np.ones(N_runs),1e5*np.ones(N_runs),320*np.ones(N_runs),1e5*np.ones(N_runs),320*np.ones(N_runs),1e5*np.ones(N_runs)))
m1 = []
m2 = []

#draw masses from the posterior
burnin = 1000
ndim = 13
filename = '../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy'
samples = np.load(filename)[:, burnin:, :].reshape((-1, ndim))
for theta in samples[np.random.randint(len(samples), size=len(K))]:
    m1.append(theta[0]/theta[10])
    m2.append(theta[1]/theta[10])

random.seed()
runs = []
path = 'output/'
for i in xrange(0,len(K)):
    seed = int(1000*random.random())
    name = path+'min%.2e_mout%.2e_MR%.1e_K%.1e_sd%d'%(m1[i],m2[i],mig_rate[i],K[i],seed)
    runs.append((m1[i],m2[i],mig_rate[i],K[i],seed,name))

os.system('make')

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3])+' '+str(pars[4])+' '+str(pars[5]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=25)
    pool.map(execute, runs)
    pool.close()
    pool.join()
