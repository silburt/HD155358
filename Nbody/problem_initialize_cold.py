#This macro is going to explore the parameter space of migration rate and K to try and find the set of parameters that best match the observed libration amplitudes of HD155358.
#From Goldreich & Schlichting (2014) we know that there is an equilibrium eccentricity associated with each tau_e/tau_n ratio. There's also a max speed the outer planet can migrate without skipping over the resonance or undergoing overstable libration.

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

N_runs = 1000

random.seed()
runs = []
mig_rate = random.sample(np.round(np.logspace(2,7,N_runs)), N_runs)
K = np.logspace(-1,4,N_runs)
path = 'output/'
for i in xrange(0,N_runs):
    seed = int(1000*random.random())
    name = path+'migrate%.1e_K%.1e_sd%d'%(mig_rate[i],K[i],seed)
    runs.append((mig_rate[i],K[i],seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=25)
    pool.map(execute, runs)
    pool.close()
    pool.join()
