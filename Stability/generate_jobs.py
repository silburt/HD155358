#This script generates the jobs for sunnyvale and stores in the jobs/ folder. To submit the jobs to sunnyvale, do 'python submit_jobs.py'

import numpy as np
import sys
import os

Njobs = 380             #sunnyvale max = 640, scinet max = 32 (times 8 per job)
logtmax = 9
sunnyvale = 1           #1=sunnyvale, 0=scinet
dir = 'jobs2/'

#delete existing jobs
#os.system('rm %s*'%dir)

#samples
filename = '../emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy'
samples = np.load(filename)[:, 1000:, :].reshape((-1, 13))

if sunnyvale == 1:  #sunnyvale
    seeds = np.random.randint(len(samples), size=Njobs)
    for i in range(Njobs):
        th = samples[seeds[i]]
        job_name = "tmax1e%d_sunnyrun%d"%(logtmax,seeds[i])
        sh_script_name = "%s%s"%(dir,job_name)
        with open(sh_script_name, 'w') as f:
            f_head = open('job_header_sunnyvale','r')
            f.write(f_head.read())
            f_head.close()
            f.write('#PBS -N %s \n'%job_name)
            f.write('# EVERYTHING ABOVE THIS COMMENT IS NECESSARY, SHOULD ONLY CHANGE nodes,ppn,walltime and my_job_name VALUES\n')
            f.write('cd $PBS_O_WORKDIR\n')      #This will be the home Stability/ directory
            f.write('./rebound output/%s %e %e %e %e %e %e %e %e %e %e %e %e >& batch.output\n'%(job_name,2*np.pi*10**logtmax,th[0],th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10]))
        f.close()

else:               #scinet - 8 jobs per batch script
    for i in range(Njobs):
        job_name = "tmax1e%d_scinetrun%d"%(logtmax,i)
        sh_script_name = "%s%s"%(dir,job_name)
        with open(sh_script_name, 'w') as f:
            f_head = open('job_header_scinet','r')
            f.write(f_head.read())
            f_head.close()
            for j in range(8):
                seed = np.random.randint(len(samples))
                th = samples[seed]
                f.write('(./rebound output/%s_sd%d %e %e %e %e %e %e %e %e %e %e %e %e) &\n'%(job_name,seed,2*np.pi*10**logtmax,th[0],th[1],th[2],th[3],th[4],th[5],th[6],th[7],th[8],th[9],th[10]))
            f.write('wait')
        f.close()
