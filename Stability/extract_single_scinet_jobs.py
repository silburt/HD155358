#This macro takes all the scinet job files and extracts them into single job files for ease of plot_stability.py.
#if repackage_singles == 1, then find the jobs that still need to be simulated longer, and repackage them into grouips of 8.
#Takes input_dir, output_dir and csv_dir as arguments

import numpy as np
import glob
import sys
import os

input_dir = sys.argv[1]     #original job files (sims grouped in 8)
output_dir = sys.argv[2]    #where you want the single job files to go

repackage_singles = 1      #take runs that need to be simulated longer and package them up

files = glob.glob('%s*'%input_dir)
for f in files:
    lines = open(f,'r').readlines()
    for l in lines:
        if './rebound' in l:
            output_name = l.split()[1].split('/')[1]
            g = open('%s%s'%(output_dir, output_name),'w')
            output = l.split(')')[0].split('(')[1]
            g.write(output)
            g.close()

#repackage jobs
if repackage_singles == 1:
    print "repackaging jobs"
    csv_dir = sys.argv[3]       #where csvs are stored
    jobs = glob.glob('%s*'%output_dir)
    fileno = 900
    njobsinfile = 0
    job_name = "tmax1e9_scinetrun%d"%(fileno)
    f = open("%s%s"%(input_dir,job_name),'w')
    f.write(open('job_header_scinet','r').read())
    for j in jobs:
        job_name = j.split(output_dir)[1]
        csv = '%s%s.csv'%(csv_dir,job_name)
        if os.path.exists(csv):
            time, dE, stable, a1, e1, l1, w1, M1, a2, e2, l2, w2, M2, phi1, phi2, phi3 = np.loadtxt(open(csv, 'r'), delimiter=',', unpack=True)
            if stable[-1] == 1 and time[-1]<6.23e9:
                j_ = open(j,'r').read()
                f.write('(%s) &\n'%j_)
                njobsinfile += 1
        if njobsinfile == 8:
            njobsinfile = 0
            f.write('wait')
            f.close()
            fileno += 1
            job_name = "tmax1e9_scinetrun%d"%(fileno)
            f = open("%s%s"%(input_dir,job_name),'w')
            f.write(open('job_header_scinet','r').read())
