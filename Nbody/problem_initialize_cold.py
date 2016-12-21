#This macro is going to explore the parameter space of migration rate and K to try and find the set of parameters that best match the observed libration amplitudes of HD155358.
#From Goldreich & Schlichting (2014) we know that there is an equilibrium eccentricity associated with each tau_e/tau_n ratio. There's also a max speed the outer planet can migrate without skipping over the resonance or undergoing overstable libration.

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

N_runs = 12

random.seed()
runs = []
mig_rate = np.round(np.logspace(3,6,N_runs))
K = 10
path = 'output/'
for i in xrange(0,N_runs):
    seed = int(1000*random.random())
    name = path+'migrate%.2e_K%d_sd%d'%(mig_rate[i],K,seed)
    runs.append((mig_rate[i],K,seed,name))

os.system('make')

length = len(runs)

def execute(pars):
    #print str(pars[0]), str(pars[1]), str(pars[2]), str(pars[3])
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2])+' '+str(pars[3]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=6)
    pool.map(execute, runs)
    pool.close()
    pool.join()
