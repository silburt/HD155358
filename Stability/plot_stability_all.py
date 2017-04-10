#This plots the stability time vs. initial conditions for HD155358.
#The original job files (i.e having the initial conditions) must go in dir/jobs/ folder
#If running scinet jobs, since those jobs files are in batches of 8, you first need to run extract_single_scinet_jobs.py

import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
from progress.bar import Bar

#arguments
dir=str(sys.argv[1])
scatter = 0

#get initial theta values from each jobs file
def get_ini(dir,f):
    job_name = f.split(dir)[1].split('.csv')[0]
    d = np.asarray(open('%sjobs/%s'%(dir,job_name), 'r').readlines()[-1].split()[3:14]).astype('float')
    d[4],d[5] = np.sqrt(d[4]**2+d[6]**2), np.sqrt(d[5]**2+d[7]**2)      #h1,h2 -> e1,e2
    d[6],d[7] = np.arccos(d[6]/d[4]), np.arccos(d[7]/d[5])              #k1,k2 -> w1,w2
    M1, M2 = d[8]-d[6], d[9]-d[7]
    d[8], d[9] = M1, M2                                                 #replace l_1 and l_2 with M_1 and M_2
    return d

def get_MAP():
    filename = '../emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy'
    samples = np.load(filename)[:, 1000:, :].reshape((-1, 13))
    d = np.percentile(samples, 50, axis=0)[:11]
    d[4],d[5] = np.sqrt(d[4]**2+d[6]**2), np.sqrt(d[5]**2+d[7]**2)      #h1,h2 -> e1,e2
    d[6],d[7] = np.arccos(d[6]/d[4]), np.arccos(d[7]/d[5])              #k1,k2 -> w1,w2
    M1, M2 = d[8]-d[6], d[9]-d[7]
    d[8], d[9] = M1, M2                                                 #replace l_1 and l_2 with M_1 and M_2
    return d

#get data
files = glob.glob('%s*.csv'%dir)
t_arr, params = [], []
bar = Bar('Processing', max=len(files))
for f in files:
    try:
        data = np.loadtxt(open(f, 'r'), delimiter=',')
        if data[-1,0] > 6.25e9 or data[-1,2] == 0:  #completed simulations
            ini = get_ini(dir,f)
            params.append(ini)
            t_arr.append(data[-1,0]/(2*np.pi))
    except:
        print "couldn't process %s"%f
    bar.next()
bar.finish()

theta = ['m_1sini','m_2sini','a_1','a_2','e_1','e_2','w_1','w_2','M_1','M_2','sini']
params = np.asarray(zip(*params))
t_arr = np.log10(np.asarray(t_arr))
MAP = get_MAP()                     #get MAP values

#plot data
size=30
fontsize=20
alpha = 0.2
fig, ax = plt.subplots(nrows=4, ncols=3, figsize=(16,15), sharey=True)
plt.subplots_adjust(wspace = 0, hspace=0.3)
ax = ax.flatten()
max = 0
for i in range(len(theta)):
    if scatter == 1:
        ax[i].scatter(params[i],t_arr, s=size, alpha=alpha, lw=0)
        ax[i].set_ylim([4,9.5])
        max = 9.6
    else:
        vals = params[i]
        x = [vals[t_arr >= 8.99], vals[t_arr < 8.99]]
        n, bins, patches = ax[i].hist(x, 20, histtype='bar', stacked=True, label=['stable','unstable'])
        nmax = np.max(n)
        max = np.max((nmax, max))
    ax[i].set_xlabel('$%s$'%theta[i], fontsize=fontsize)
    ax[i].plot([MAP[i], MAP[i]], [0,max], 'k--', lw=3)
    xticks = ax[i].xaxis.get_major_ticks()
    xticks[-1].label1.set_visible(False)
    #xticks[0].label1.set_visible(False)

#final figure fixes
fig.delaxes(ax[11])
ax[10].set_xlim([0,1])
if scatter == 0:
    ax[0].legend(loc='upper left', fontsize=11)

#save
ext = 'histogram'
if scatter == 1:
    ext = 'scatter'
plt.savefig('%sstability_plot_%s.pdf'%(dir,ext))

#Percent that are stable
N_samples = float(len(t_arr))
stable = float(len(t_arr[t_arr>8.99]))/N_samples
print "percent of stable systems over 1e9 years = %f +/- %f"%(stable,np.sqrt(N_samples)/N_samples)
print "Total samples used: %d"%len(t_arr)
