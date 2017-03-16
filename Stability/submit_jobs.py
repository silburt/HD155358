import os
import glob

os.system('make clean && make')

jobs_dir = 'jobs/'
files = glob.glob('%s/*'%jobs_dir)
for f in files:
    job_name = f.split(jobs_dir)[1]
    os.system('mv %s %s'%(f, job_name))
    os.system('qsub %s'%job_name)           #submit job in home directory
    os.system('mv %s %s'%(job_name,f))
