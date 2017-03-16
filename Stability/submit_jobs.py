import os
import glob

jobs_dir = 'jobs/'

os.system('make clean && make')
os.system('mv rebound %s.'%jobs_dir)
os.system('cd %s'%jobs_dir)

files = glob.glob()
for f in files:
    #job_name = f.split(jobs_dir)[1]
    #os.system('mv %s %s'%(f, job_name))
    os.system('qsub %s'%f)           #submit job in home directory
    #os.system('mv %s %s'%(job_name,f))
