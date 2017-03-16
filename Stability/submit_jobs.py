import os
import glob

files = glob.glob('jobs/*')
for f in files:
    os.system('qsub %s'%f)
