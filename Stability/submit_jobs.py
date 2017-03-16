import os
import glob

jobs_dir = 'jobs/'

os.system('make clean && make')
#os.system('mv rebound %s.'%jobs_dir)
#os.chdir(jobs_dir)

files = glob.glob('%s/*'%jobs_dir)
for f in files:
    job_name = f.split(jobs_dir)[1]
    os.system('mv %s %s'%(f, job_name))
    os.system('qsub %s'%job_name)
    os.system('mv %s %s'%(job_name,f))
