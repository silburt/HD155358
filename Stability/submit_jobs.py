#This script submits jobs to the sunnyvale cluster (and resubmits jobs that stopped because of the 48hour walltime limit).

import os
import os.path
import glob
import numpy as np

def submit_job(f, job_name):
    os.system('mv %s %s'%(f, job_name))
    os.system('qsub %s'%job_name)
    os.system('mv %s %s'%(job_name,f))

jobs_dir = 'jobs/'

os.system('make clean && make')

files = glob.glob('%s/*'%jobs_dir)
Njobs_counter = 0
for f in files:
    job_name = f.split(jobs_dir)[1]
    
    #if integration has already been run, check if was still stable
    csv = 'output/%s.csv'%job_name
    if os.path.exists(csv):
        time, dE, stable, a1, e1, l1, w1, M1, a2, e2, l2, w2, M2, phi1, phi2, phi3 = np.loadtxt(open(csv, 'r'), delimiter=',', unpack=True)
        if stable[-1] == 1 and time[-1]<6.26e9:
            print 'restarting job: %s'%job_name
            submit_job(f, job_name)
            Njobs_counter += 1
    else:
        submit_job(f, job_name)     #submitting job for the first time
        Njobs_counter += 1

    if Njobs_counter >= 635:        #640 single-node job limit for sunnyvale
        break

print 'submitted %d jobs'%Njobs_counter
