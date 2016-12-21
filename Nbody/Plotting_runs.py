#Testing how elapsed time scales with N_planetesimals

import multiprocessing as mp
import os
import sys
import time
import random
import numpy as np

#Specify what runs you want *****************************
N_runs = 12

random.seed()
runs = []
N_planetesimals = np.round(np.logspace(2.7,4.905,N_runs))
name = 'output/Aug24'
for i in xrange(0,N_runs):
    Np = "{:.0f}".format(N_planetesimals[i])
    seed = "{:.0f}".format(int(1000*random.random()))
    runs.append((name,Np,seed))

os.system('make')

length = len(runs)

def execute(pars):
    os.system('./rebound '+str(pars[0])+' '+str(pars[1])+' '+str(pars[2]))

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=length)
    pool.map(execute, runs)
    pool.close()
    pool.join()
