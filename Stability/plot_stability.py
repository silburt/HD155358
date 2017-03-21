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

#get initial theta values from jobs file
def get_ini(dir,f,convert):
    job_name = f.split(dir)[1].split('.csv')[0]
    d = np.asarray(open('%sjobs/%s'%(dir,job_name), 'r').readlines()[-1].split()[3:]).astype('float')
    d[0] /= d[-1] #m1sini -> m1
    d[1] /= d[-1] #m2sini -> m2
    if convert == 1:
        d[4],d[5] = np.sqrt(d[4]**2+d[6]**2), np.sqrt(d[5]**2+d[7]**2)      #h1,h2 -> e1,e2
        d[6],d[7] = np.arccos(d[6]/d[4]), np.arccos(d[7]/d[5])              #k1,k2 -> w1,w2
        d[8],d[9] = d[8]-d[6], d[9]-d[7]                                    #l1,l2 -> M1,M2
    return d

#dictionary mapping theta to index values
theta = ['m_1','m_2','a_1','a_2','e_1','e_2','w_1','w_2','M_1','M_2','sini']
orb = dict()
for index, orb_ele in enumerate(theta):
    orb[orb_ele] = index

#get data
corb = 0    #convert to orbital elements from Pal
if orb[x1] > 3 or orb[y1] > 3:
    corb = 1
x1_arr, y1_arr, x2_arr, y2_arr, t_arr = [], [], [], [], []
bar = Bar('Processing', max=len(files[:100]))
for f in files[:100]:
    #time, dE, stable, a1, e1, l1, w1, M1, a2, e2, l2, w2, M2, phi1, phi2, phi3 = np.loadtxt(open(f, 'r'), delimiter=',', unpack=True)
    data = np.loadtxt(open(f, 'r'), delimiter=',')
    ini = get_ini(dir,f,corb)
    x1_arr.append(ini[orb[x1]]), y1_arr.append(ini[orb[y1]]), x2_arr.append(ini[orb[x2]]), y2_arr.append(ini[orb[y2]]), t_arr.append(data[-1,0]/(2*np.pi))
    bar.next()
bar.finish()

#plot data
size=15
colorbar = 'winter'
fontsize=16
fig, ax = plt.subplots(nrows=2, ncols=1, figsize=(8,10))

ax[0].scatter(x1_arr,y1_arr,c=np.log10(np.asarray(t_arr)), s=size, cmap=colorbar, lw=0, label='simulated points')
ax[0].set_xlabel(r'$%s$'%x1)
ax[0].set_ylabel(r'$%s$'%y1)
sc = ax[1].scatter(x2_arr,y2_arr,c=np.log10(np.asarray(t_arr)), s=size, cmap=colorbar, lw=0, label='simulated points')
ax[1].set_xlabel(r'$%s$'%x2)
ax[1].set_ylabel(r'$%s$'%y2)

fig.subplots_adjust(right=0.8)
cbar_ax = fig.add_axes([0.85, 0.15, 0.05, 0.7])
fig.colorbar(sc, label=r'log10($t_{stable}$)',cax=cbar_ax)
plt.savefig('%sstability_plot.pdf'%dir)
