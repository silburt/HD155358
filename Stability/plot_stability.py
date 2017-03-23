#This plots the stability time vs. initial conditions for HD155358. These use runs from Sunnyvale

import numpy as np
import glob
import sys
import matplotlib.pyplot as plt
from progress.bar import Bar

#initial conditions
dir=str(sys.argv[1])
files = glob.glob('%s*.csv'%dir)
x1,y1 = 'a_1','e_1'
x2,y2 = 'a_2','e_2'
x3,y3 = 'a_2','w_2'
x4,y4 = 'm_1','m_2'

#get initial theta values from jobs file
def get_ini(dir,f):
    job_name = f.split(dir)[1].split('.csv')[0]
    d = np.asarray(open('%sjobs/%s'%(dir,job_name), 'r').readlines()[-1].split()[3:14]).astype('float')
    d[0] /= d[10]                                                       #m1sini -> m1
    d[1] /= d[10]                                                       #m2sini -> m2
    d[4],d[5] = np.sqrt(d[4]**2+d[6]**2), np.sqrt(d[5]**2+d[7]**2)      #h1,h2 -> e1,e2
    d[6],d[7] = np.arccos(d[6]/d[4]), np.arccos(d[7]/d[5])              #k1,k2 -> w1,w2
    d = np.append(d, (d[8]-d[6],d[9]-d[7]))                             #add M1,M2 (keep l1, l2)
    return d

def get_MAP():
    filename = '../emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy'
    samples = np.load(filename)[:, 1000:, :].reshape((-1, 13))
    d = np.percentile(samples, 50, axis=0)[:11]
    d[4],d[5] = np.sqrt(d[4]**2+d[6]**2), np.sqrt(d[5]**2+d[7]**2)      #h1,h2 -> e1,e2
    d[6],d[7] = np.arccos(d[6]/d[4]), np.arccos(d[7]/d[5])              #k1,k2 -> w1,w2
    d = np.append(d, (d[8]-d[6],d[9]-d[7]))                             #add M1,M2 (keep l1, l2)
    return d

#dictionary mapping theta to index values
theta = ['m_1','m_2','a_1','a_2','e_1','e_2','w_1','w_2','l_1','l_2','sini','M_1','M_2']
orb = dict()
for index, orb_ele in enumerate(theta):
    orb[orb_ele] = index

#get data
t_arr, x1_arr, y1_arr, x2_arr, y2_arr, x3_arr, y3_arr, x4_arr, y4_arr = [], [], [], [], [], [], [], [], []
bar = Bar('Processing', max=len(files))
for f in files:
    data = np.loadtxt(open(f, 'r'), delimiter=',')
    ini = get_ini(dir,f)
    x1_arr.append(ini[orb[x1]]), y1_arr.append(ini[orb[y1]]), x2_arr.append(ini[orb[x2]]), y2_arr.append(ini[orb[y2]])
    x3_arr.append(ini[orb[x3]]), y3_arr.append(ini[orb[y3]]), x4_arr.append(ini[orb[x4]]), y4_arr.append(ini[orb[y4]])
    t_arr.append(data[-1,0]/(2*np.pi))
    bar.next()
bar.finish()

#plot data
size=15
colorbar = 'winter'
fontsize=20
fig, ax = plt.subplots(nrows=2, ncols=2, figsize=(16,15))
ax = ax.flatten()
ax[0].scatter(x1_arr,y1_arr,c=np.log10(np.asarray(t_arr)), s=size, cmap=colorbar, lw=0)
ax[0].set_xlabel(r'$%s$'%x1,fontsize=fontsize)
ax[0].set_ylabel(r'$%s$'%y1,fontsize=fontsize)
ax[1].scatter(x2_arr,y2_arr,c=np.log10(np.asarray(t_arr)), s=size, cmap=colorbar, lw=0)
ax[1].set_xlabel(r'$%s$'%x2,fontsize=fontsize)
ax[1].set_ylabel(r'$%s$'%y2,fontsize=fontsize)
ax[2].scatter(x3_arr,y3_arr,c=np.log10(np.asarray(t_arr)), s=size, cmap=colorbar, lw=0)
ax[2].set_xlabel(r'$%s$'%x3,fontsize=fontsize)
ax[2].set_ylabel(r'$%s$'%y3,fontsize=fontsize)
sc = ax[3].scatter(x4_arr,y4_arr,c=np.log10(np.asarray(t_arr)), s=size, cmap=colorbar, lw=0)
ax[3].set_xlabel(r'$%s$'%x4,fontsize=fontsize)
ax[3].set_ylabel(r'$%s$'%y4,fontsize=fontsize)

#overplot MAP values
MAP = get_MAP()
ax[0].scatter(MAP[orb[x1]],MAP[orb[y1]],marker='+',color='black',s=1000)
ax[1].scatter(MAP[orb[x2]],MAP[orb[y2]],marker='+',color='black',s=1000)
ax[2].scatter(MAP[orb[x3]],MAP[orb[y3]],marker='+',color='black',s=1000)
ax[3].scatter(MAP[orb[x4]],MAP[orb[y4]],marker='+',color='black',s=1000)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
cb = fig.colorbar(sc, label=r'log10($t_{stable}$)',cax=cbar_ax)
cb.ax.set_ylabel(cb.ax.get_ylabel(), fontsize=fontsize)
plt.savefig('%sstability_plot.pdf'%dir)