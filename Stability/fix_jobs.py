import glob
import os

jobs = glob.glob('jobs/*')
for f in jobs:
    fileOpen = open(f, 'a')
    lines = open(f, 'r').readlines()
    new_last_line = (lines[-1].rstrip() + " >& batch.output")
    lines[-1] = new_last_line
    open(f, 'w').writelines(lines)
