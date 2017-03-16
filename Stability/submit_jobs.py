import os
import glob

jobs_dir = 'jobs/'

os.system('make clean && make')

files = glob.glob('%s/*'%jobs_dir)
for f in files:
    job_name = f.split(jobs_dir)[1]
    os.system('mv %s %s'%(f, job_name))
    os.system('qsub %s'%job_name)
    os.system('mv %s %s'%(job_name,f))
