import os
import glob

jobs_dir = 'jobs/'

os.system('make clean && make')
os.system('mv rebound %s.'%jobs_dir)
os.chdir(jobs_dir)

files = glob.glob('tmax*')
for f in files:
    os.system('qsub %s'%f)

