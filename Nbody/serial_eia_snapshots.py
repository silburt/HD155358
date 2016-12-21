#This macro shows the evolution of semimajor axis, eccentricity and inclination for different snapshots in time. Done in a serial fashion.

import sys
import matplotlib.pyplot as plt
import numpy as np
import glob
import re

def natural_key(string_):
    return [int(s) if s.isdigit() else s for s in re.split(r'(\d+)', string_)]

dir = sys.argv[1]
files = glob.glob(dir+'*_eiasnapshot_t=0.txt')

print 'getting '+str(len(files))+' snapshots'
for f in files:
    f = f.split('0.txt')[0]
    snapshots = glob.glob(f+'*')
    snapshots = sorted(snapshots, key = natural_key)
    Nsnapshots = len(snapshots)

    snapshot_name = re.sub('_t=0\.txt$', '', snapshots[0])
    if Nsnapshots > 1:
        fig, axes = plt.subplots(nrows=Nsnapshots, ncols=1, figsize=(12,12), sharex=True)
        for j in xrange(0,Nsnapshots):
            fos = open(snapshots[j], 'r')
            time, id, a, e, inc, rdist, m = np.loadtxt(fos, delimiter=',',unpack=True)
            axes[j].plot(a,e,'.')
            axes[j].plot(a[0],e[0],'o',color='red')
            axes[j].plot(a[1],e[1],'o',color='orange')
            axes[j].set_title('t = '+str(time[0])+' years, N='+str(len(a)))
            axes[j].yaxis.set_ticks(np.arange(0,max(e),0.1))
        axes[Nsnapshots-1].set_ylabel('$e$')
        axes[Nsnapshots-1].set_xlabel('(a (AU)')
        plt.savefig(snapshot_name+'s.png')
