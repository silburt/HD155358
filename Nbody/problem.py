#This macro is going to explore the parameter space of migration rate and K to try and find the set of parameters that best match the observed libration amplitudes of HD155358.
#From Goldreich & Schlichting (2014) we know that there is an equilibrium eccentricity associated with each tau_e/tau_n ratio. There's also a max speed the outer planet can migrate without skipping over the resonance or undergoing overstable libration.

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

N_runs = 1

#draw masses from the posterior
m1 = []
m2 = []
sini = []
burnin = 1000
ndim = 13
filename = '../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy'
samples = np.load(filename)[:, burnin:, :].reshape((-1, ndim))
for theta in samples[np.random.randint(len(samples), size=N_runs)]:
    m1.append(theta[0])
    m2.append(theta[1])
    sini.append(theta[10])

random.seed()
runs = []
mig_rate = random.sample(np.round(np.logspace(2,7,N_runs)), N_runs)
K1 = random.sample(np.logspace(-1,4,N_runs), N_runs)
K2 = np.logspace(-1,4,N_runs)
path = 'output/'
for i in xrange(0,N_runs):
    seed = int(1000*random.random())
    name = path+'taueinner_migrate%.1e_Kin%.1e_Kout%.1e_sd%d'%(mig_rate[i],K1[i],K2[i],seed)
    runs.append((m1[i],m2[i],sini[i],mig_rate[i],K1[i],K2[i],seed,name))

os.system('make')

def execute(pars):
    os.system('./rebound %f %f %f %f %f %f %d %s'%pars)

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=25)
    pool.map(execute, runs)
    pool.close()
    pool.join()