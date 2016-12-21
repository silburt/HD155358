#Compute all the orbits.py for a given folder

import multiprocessing as mp
import os
import sys
import time
import random
import glob


#Specify directory
dir = sys.argv[1]
orbital_param = 'P'
N_cutoffs = [6310,12589,25119,50119,100000]
params = []
for N in N_cutoffs:
    params.append((dir,str(N),orbital_param))
N_runs = len(params)

def execute(p):
    os.system('python NpvEvolution.py '+p[0]+' '+p[1]+' '+p[2])
    print 'finished process'

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=N_runs)
    print 'computing '+str(N_runs)+' NpvEvolution.py analyses'
    pool.map(execute, params)
    pool.close()
    pool.join()