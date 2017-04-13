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

def fixr(a):    #fix range to 0-2pi
    while a >= 2*np.pi:
        a -= 2*np.pi
    while a < 0:
        a += 2*np.pi
    return a

def get_theta(theta,dir='',f=''):
    m1sini,m2sini,a1,a2,h1,h2,k1,k2,l1,l2,sini,J2,off = theta
    e1, e2 = np.sqrt(h1**2+k1**2), np.sqrt(h2**2+k2**2)
    w1, w2 = fixr(np.arctan2(h1,k1)), fixr(np.arctan2(h2,k2))
    M1, M2 = fixr(l1-k1), fixr(l2-k2)

    #Check that theta parameters (from MCMC chain) match the csv (stability sim) file
    if dir != '':
        job_name = f.split(dir)[1].split('.csv')[0]
        d = np.asarray(open('%sjobs/%s'%(dir,job_name), 'r').readlines()[-1].split()[3:14]).astype('float')
        if np.sum(theta[0:11]-d) > 1e-5:
            print "parameter mismatch for %s"%f
    return [m1sini,m2sini,sini,a1,a2,off,e1,e2,J2**0.5,w1,w2,0,M1,M2]

#get data
samples = np.load('../emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy')[:, 1000:, :].reshape((-1, 13))
files = glob.glob('%s*.csv'%dir)
t_arr, params = [], []
bar = Bar('Processing', max=len(files))
for f in files:
    try:
        data = np.loadtxt(open(f, 'r'), delimiter=',')
        if np.any(np.isnan(data)):
            print 'NaN in %s'%f
        if data[-1,0] > 6.25e9 or data[-1,2] == 0:  #completed simulations
            try:
                seed = int(f.split(dir)[1].split('.csv')[0].split('sd')[1])
            except:
                seed = int(f.split(dir)[1].split('.csv')[0].split('sunnyrun')[1])
            ini = get_theta(samples[seed],dir,f)
            params.append(ini)
            t_arr.append(data[-1,0]/(2*np.pi))
    except:
        print "couldn't process %s"%f
    bar.next()
bar.finish()

theta = ['m_1\mathrm{sin}i \ (\mathrm{M_J})','m_2\mathrm{sin}i \ (\mathrm{M_J})','sini','a_1 \ \mathrm{(AU)}','a_2 \ \mathrm{(AU)}','\gamma \ \mathrm{(m/s)}','e_1','e_2','J \ \mathrm{(m/s)}','w_1','w_2','junk','M_1','M_2']
params = np.asarray(zip(*params))
t_arr = np.log10(np.asarray(t_arr))
MAP = get_theta(np.percentile(samples,50,axis=0))

#plot data
size=30
fontsize=20
alpha = 0.2
fig, ax = plt.subplots(nrows=5, ncols=3, figsize=(16,15), sharey=True)
plt.subplots_adjust(wspace = 0, hspace=0.4)
ax = ax.flatten()
max = 0
for i in range(len(theta)):
    if i==11:   #skip 11th plot, to be deleted later
        continue
    vals = params[i]
    x = [vals[t_arr >= 8.99], vals[t_arr < 8.99]]   #split into stable/unstable systems
    n, bins, patches = ax[i].hist(x, 20, histtype='bar', stacked=True, color=['lime','red'], label=['stable','unstable'], edgecolor='none')
    ax[i].set_xlabel(r'$%s$'%theta[i], fontsize=fontsize)
    xticks = ax[i].xaxis.get_major_ticks()
    xticks[-1].label1.set_visible(False)
    #xticks[0].label1.set_visible(False)

ymin, ymax = ax[i].get_ylim()
for i in range(len(theta)):
    ax[i].plot([MAP[i], MAP[i]], [0,ymax], 'k--', lw=3)

#final figure fixes
fig.delaxes(ax[11])
fig.delaxes(ax[14])
ax[2].set_xlim([0,1])   #sini
ax[0].legend(loc='upper left', fontsize=11)

#save
plt.savefig('%sstability_plot_histogram.pdf'%dir)

#Percent that are stable
N_samples = float(len(t_arr))
stable = float(len(t_arr[t_arr>8.99]))/N_samples
print "percent of stable systems over 1e9 years = %f +/- %f"%(stable,np.sqrt(N_samples)/N_samples)
print "Total samples used: %d"%len(t_arr)


####################
######old code######
####################
#get initial theta values from each jobs file
#def get_ini(dir,f):
#    job_name = f.split(dir)[1].split('.csv')[0]
#    d = np.asarray(open('%sjobs/%s'%(dir,job_name), 'r').readlines()[-1].split()[3:14]).astype('float')
#    w1, w2 = np.arctan2(d[4],d[6]), np.arctan2(d[5],d[7])           #k1,k2 -> w1,w2
#    M1, M2 = d[8]-d[6], d[9]-d[7]
#    d[4],d[5] = np.sqrt(d[4]**2+d[6]**2), np.sqrt(d[5]**2+d[7]**2)  #h1,h2 -> e1,e2
#    d[6],d[7],d[8],d[9] = fixr(w1),fixr(w2),fixr(M1),fixr(M2)       #replace k1,k2->w1,w2 and l_1,l_2->M_1,M_2
#    return d

#def get_MAP():
#    filename = '../emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy'
#    samples = np.load(filename)[:, 1000:, :].reshape((-1, 13))
#    d = np.percentile(samples, 50, axis=0)[:11]
#    w1, w2 = np.arctan2(d[4],d[6]), np.arctan2(d[5],d[7])           #k1,k2 -> w1,w2
#    M1, M2 = d[8]-d[6], d[9]-d[7]
#    d[4],d[5] = np.sqrt(d[4]**2+d[6]**2), np.sqrt(d[5]**2+d[7]**2)  #h1,h2 -> e1,e2
#    d[6],d[7],d[8],d[9] = fixr(w1),fixr(w2),fixr(M1),fixr(M2)       #replace k1,k2->w1,w2 and l_1,l_2->M_1,M_2
#    return d
