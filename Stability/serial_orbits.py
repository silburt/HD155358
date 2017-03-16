#Compute all the orbits.py for a given folder

import multiprocessing as mp
import os
import sys
import time
import random
import glob

#ext = '_Mpl10mp'    #what goes right before the .txt

#Specify directory
dir = sys.argv[1]
files = glob.glob(dir+'*.csv')
N = len(files)

def execute(filename):
    os.system('python orbits.py '+filename)
    print 'finished process'

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    pool = mp.Pool(processes=10)
    print 'computing '+str(N)+' orbit.py analyses'
    pool.map(execute, files)
    pool.close()
    pool.join()
    #os.system('mv %s/*_orbit.png %s/orbits/.'%(dir, dir))
