#This macro plots the results in a 2d heatmap: K vs. migration rate, with the color being the various amplitudes (period ratio, eccentricity, etc).

import numpy as np
import pandas as pd
import glob
import sys
import matplotlib.pyplot as plt
from progress.bar import Bar
import rebound
from mpl_toolkits.mplot3d import Axes3D

def check_resonance(theta, times):
    m1sini,m2sini,aa1,aa2,h1,h2,k1,k2,lambda1,lambda2,sini = theta
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
    mJ = 9.543e-4                        #Jupiter mass -> solar mass
    a1, e1, a2, e2, Pr, phi1, phi2, phi3 = [], [], [], [], [], np.empty(0), np.empty(0), np.empty(0)
    
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.add(m=0.92)                      #add the star
    sim.dt = 2*np.pi* aa1**(1.5) / 25.    #dt = 25 steps per orb per of inner planet
    sim.add(m=m1sini*mJ/sini,a=aa1,l=lambda1,h=h1,k=k1)
    sim.add(m=m2sini*mJ/sini,a=aa2,l=lambda2,h=h2,k=k2)
    sim.move_to_com()
    for t in times:
        sim.integrate(t)
        o1 = sim.particles[1].calculate_orbit(sim.particles[0])
        o2 = sim.particles[2].calculate_orbit(sim.particles[0])
        a1.append(o1.a),e1.append(o1.e),a2.append(o2.a),e2.append(o2.e), Pr.append(o2.P/o1.P)
        phi1 = np.append(phi1,2.*o2.l - o1.l - o1.pomega)
        phi2 = np.append(phi1,2.*o2.l - o1.l - o2.pomega)
        phi3 = np.append(phi1,o1.pomega - o2.pomega)
    while np.any(phi1>2*np.pi) or np.any(phi2>2*np.pi) or np.any(phi3>2*np.pi):
        phi1[phi1>2*np.pi] -= 2*np.pi
        phi2[phi2>2*np.pi] -= 2*np.pi
        phi3[phi3>2*np.pi] -= 2*np.pi
    while np.any(phi1<0) or np.any(phi2<0) or np.any(phi3<0):
        phi1[phi1<0] += 2*np.pi
        phi2[phi2<0] += 2*np.pi
        phi3[phi3<0] += 2*np.pi
    return a1, e1, a2, e2, Pr, phi1, phi2, phi3

def get_phi_amplitude(phi, endpoints):
    A1 = np.max(phi[-endpoints:]) - np.min(phi[-endpoints:])  #by default, these angles are between 0-2pi
    phi[phi > np.pi] -= 2*np.pi
    A2 = np.max(phi[-endpoints:]) - np.min(phi[-endpoints:])
    return min([A1,A2])

#Specify directory
dir = sys.argv[1]
files = glob.glob(dir+'*.txt')
N = len(files)
i=0
while i < N:            #just want the main .txt files
    f = files[i]
    string = f.split("_")
    if string[-1]=="info.txt" or string[-1]=="elapsedtime.txt" or string[-2]=="eiasnapshot":
        files.remove(files[i])
        N -= 1
    else:
        i += 1

bar = Bar('Processing', max=N)
endpoints = 400
A_P = []    #period amplitude
A_e1 = []   #eccentricity amplitude of planet 1
A_e2 = []   #eccentricity amplitude of planet 2
A_a1 = []   #semi-major axis amplitude of planet 1
A_a2 = []   #semi-major axis amplitude of planet 2
MR = []     #migration rate
K1 = []     #inner planet damping constant Lee & Peale (2002)
K2 = []     #outer planet damping constant Lee & Peale (2002)
A_phi1, A_phi2, A_phi3 = [], [], []
for i,f in enumerate(files):
    fos = open(f, 'r')
    #time, dE, N, mig_rate, dampratio, a1, e1, a2, e2, phi1, phi2, phi3 = np.loadtxt(fos, delimiter=',', unpack=True)
    time, dE, N, mig_rate, damp1, damp2, a1, e1, a2, e2, phi1, phi2, phi3, m1, m2, taua1, taue1, taua2, taue2 = np.loadtxt(fos, delimiter=',', unpack=True)
    if np.any(e2>1) or np.any(e2<0) or np.any(e1>1) or np.any(e1<0) or np.any(np.abs((a2/a1)**1.5 - 2)>0.3):
        pass
        #print "\nerror in %s, skipping file."%f
    else:
        Pratio = (a2[-endpoints:]/a1[-endpoints:])**1.5
        A_P.append((np.max(Pratio) - np.min(Pratio))/2.)
        A_e1.append((np.max(e1[-endpoints:]) - np.min(e1[-endpoints:]))/2.)
        A_e2.append((np.max(e2[-endpoints:]) - np.min(e2[-endpoints:]))/2.)
        A_a1.append((np.max(a1[-endpoints:]) - np.min(a1[-endpoints:]))/2.)
        A_a2.append((np.max(a2[-endpoints:]) - np.min(a2[-endpoints:]))/2.)
        A_phi1.append(get_phi_amplitude(phi1, endpoints))
        A_phi2.append(get_phi_amplitude(phi2, endpoints))
        A_phi3.append(get_phi_amplitude(phi3, endpoints))
        MR.append(mig_rate[0])
        K1.append(damp1[0])
        K2.append(damp2[0])
    bar.next()
bar.finish()

cols = ["P", "e1", "e2", "a1", "a2", "phi1", "phi2", "phi3"]
d1 = pd.DataFrame(zip(A_P, A_e1, A_e2, A_a1, A_a2, A_phi1, A_phi2, A_phi3), columns=cols)

print "Finished analyzing simulated runs. Now analyzing samples from posterior to plot"
#draw samples from emcee chain
Ndraws = 100
tmax = 2000
Npts = 1000
times = np.linspace(0,tmax,Npts)

burnin = 1000
ndim = 13
filename = '../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy'
samples = np.load(filename)[:, burnin:, :].reshape((-1, ndim))
bar = Bar('Processing', max=Ndraws)
A_Ps, A_e1s, A_e2s, A_a1s, A_a2s, A_phi1s, A_phi2s, A_phi3s = [], [], [], [], [], [], [], []
endpoints = 250
for theta in samples[np.random.randint(len(samples), size=Ndraws)]:
    a1, e1, a2, e2, Pratio, phi1, phi2, phi3 = check_resonance(theta[:11],times)
    A_Ps.append((np.max(Pratio) - np.min(Pratio))/2.)
    A_e1s.append((np.max(e1[-endpoints:]) - np.min(e1[-endpoints:]))/2.)
    A_e2s.append((np.max(e2[-endpoints:]) - np.min(e2[-endpoints:]))/2.)
    A_a1s.append((np.max(a1[-endpoints:]) - np.min(a1[-endpoints:]))/2.)
    A_a2s.append((np.max(a2[-endpoints:]) - np.min(a2[-endpoints:]))/2.)
    A_phi1s.append(get_phi_amplitude(phi1, endpoints))
    A_phi2s.append(get_phi_amplitude(phi2, endpoints))
    A_phi3s.append(get_phi_amplitude(phi3, endpoints))
    bar.next()
bar.finish()

d2 = pd.DataFrame(zip(A_Ps, A_e1s, A_e2s, A_a1s, A_a2s, A_phi1s, A_phi2s, A_phi3s), columns=cols)

size=25
colorbar = 'autumn'
fig = plt.figure(figsize=(8,10))
plot_x, plot_y = "phi1", "phi2"     #the x and y values to be plotted
MCMC_alpha = 0.5

#plot 1
axes = fig.add_subplot(3, 1, 1)#, projection='3d')
sc = axes.scatter(d1[plot_x],d1[plot_y],c=np.log10(K1), s=size, cmap=colorbar, lw=0, label='simulated points')
axes.scatter(d2[plot_x],d2[plot_y], color='black', s=size, cmap=colorbar, alpha=MCMC_alpha, lw=0, label='MCMC samples')
plt.colorbar(sc, ax=axes, label=r'log10($K_1$), $K_1=\tau_{a_2}/\tau_{e_1}$')
axes.legend(loc='upper left', numpoints=1, fontsize=8)
axes.set_xlabel('Amplitude %s'%plot_x)
axes.set_ylabel('Amplitude %s'%plot_y)
#axes.set_xscale('log')
#axes.set_xlim([1e-4,1])

#plot 2
axes = fig.add_subplot(3, 1, 2)
sc = axes.scatter(d1[plot_x],d1[plot_y],c=np.log10(K2), s=size, cmap=colorbar, lw=0, label='simulated points')
axes.scatter(d2[plot_x],d2[plot_y], color='black', s=size, cmap=colorbar, alpha=MCMC_alpha, lw=0, label='MCMC samples')
plt.colorbar(sc, ax=axes, label=r'log10($K_2$), $K_2=\tau_{a_2}/\tau_{e_2}$')
axes.set_xlabel('Amplitude %s'%plot_x)
axes.set_ylabel('Amplitude %s'%plot_y)

#plot 3
axes = fig.add_subplot(3, 1, 3)
sc = axes.scatter(d1[plot_x],d1[plot_y],c=np.log10(MR), s=size, cmap=colorbar, lw=0, label='simulated points')
axes.scatter(d2[plot_x],d2[plot_y], color='black', s=size, cmap=colorbar, alpha=MCMC_alpha, lw=0, label='MCMC samples')
plt.colorbar(sc, ax=axes, label='log10(Migration Rate)')
axes.set_xlabel('Amplitude %s'%plot_x)
axes.set_ylabel('Amplitude %s'%plot_y)

#save & plot
plt.savefig(dir+'amplitude_%s_%s.png'%(plot_x, plot_y))
plt.show()

'''
#3D plots
#sc = axes.scatter(A_phi1,A_phi2,A_phi3,c=np.log10(K), s=size, cmap=colorbar, lw=0, label='simulated points')
#axes.scatter(A_phi1s,A_phi2s,A_phi3s, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
#axes.set_zlabel('Amplitude Phi3')
#axes.view_init(elev = 12, azim=-91)
#sc = axes.scatter(A_phi1,A_phi2,A_phi3,c=np.log10(MR), s=size, cmap=colorbar, lw=0, label='simulated points')
#axes.scatter(A_phi1s,A_phi2s,A_phi3s, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
#axes.set_zlabel('Amplitude Phi3')
#axes.view_init(elev = 12, azim=-91)
'''

'''
#figures
size=25
colorbar = 'autumn'
fig, axes = plt.subplots(nrows=4, ncols=1, figsize=(10,10))
plt.subplots_adjust(hspace = 0.45)

sc = axes[0].scatter(A_e1, A_e2, c=np.log10(K), s=size, cmap=colorbar, lw=0, label='simulated points')
sc2 = axes[1].scatter(A_e1, A_e2, c=np.log10(MR), s=size, cmap=colorbar, lw=0, label='simulated points')
#sc3 = axes[2].scatter(A_a1, A_P, c=np.log10(K), s=size, cmap=colorbar, lw=0, label='simulated points')
#sc4 = axes[3].scatter(A_a1, A_P, c=np.log10(K), s=size, cmap=colorbar, lw=0, label='simulated points')
sc3 = axes[2].scatter(A_a1, A_a2, c=np.log10(K), s=size, cmap=colorbar, lw=0, label='simulated points')
sc4 = axes[3].scatter(A_a1, A_a2, c=np.log10(MR), s=size, cmap=colorbar, lw=0, label='simulated points')
axes[0].scatter(A_e1s, A_e2s, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
axes[1].scatter(A_e1s, A_e2s, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
#axes[2].scatter(A_a1s, A_Ps, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
#axes[3].scatter(A_a1s, A_Ps, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
axes[2].scatter(A_a1s, A_a2s, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
axes[3].scatter(A_a1s, A_a2s, color='black', s=size, cmap=colorbar, lw=0, label='MCMC samples')
axes[0].set_xlabel('e1 amplitude')
axes[0].set_ylabel('e2 amplitude')
axes[0].set_xlim([1e-3,0.4])
axes[0].set_ylim([1e-3,0.4])
axes[0].set_xscale('log')
axes[0].set_yscale('log')
axes[0].legend(loc='upper left', numpoints=1, fontsize=8)
axes[1].set_xlabel('e1 amplitude')
axes[1].set_ylabel('e2 amplitude')
axes[1].set_xlim([1e-3,0.4])
axes[1].set_ylim([1e-3,0.4])
axes[1].set_xscale('log')
axes[1].set_yscale('log')
axes[2].set_xlabel('a1 amplitude')
axes[2].set_ylabel('a2 amplitude')
axes[2].set_xlim([1e-4,0.2])
axes[2].set_ylim([1e-4,0.2])
axes[2].set_xscale('log')
axes[2].set_yscale('log')
axes[2].legend(loc='upper left', numpoints=1, fontsize=8)
axes[3].set_xlabel('a1 amplitude')
axes[3].set_ylabel('a2 amplitude')
axes[3].set_xlim([1e-4,0.2])
axes[3].set_ylim([1e-4,0.2])
axes[3].set_xscale('log')
axes[3].set_yscale('log')
plt.colorbar(sc, ax=axes[0], label=r'log10($K$), $K=\tau_e/\tau_a$')
plt.colorbar(sc2, ax=axes[1], label='log10(Migration Rate)')
plt.colorbar(sc3, ax=axes[2], label=r'log10($K$), $K=\tau_e/\tau_a$')
plt.colorbar(sc4, ax=axes[3], label='log10(Migration Rate)')
axes[0].set_title('Variable Amplitudes: 0.5*(V$_{max}$ - V$_{min}$)')
'''
'''
sc = axes[0].scatter(K, MR, c=A_P, s=size, cmap=colorbar, lw=0)
sc2 = axes[1].scatter(K, MR, c=A_e1, s=size, cmap=colorbar, lw=0)
sc3 = axes[2].scatter(K, MR, c=A_e2, s=size, cmap=colorbar, lw=0)
axes[0].set_xscale('log')
axes[0].set_xlim([1e-2,1e6])
axes[0].set_yscale('log')
axes[1].set_yscale('log')
axes[2].set_yscale('log')
axes[0].set_ylabel('Outer Planet Migration Speed')
axes[1].set_ylabel('Outer Planet Migration Speed')
axes[2].set_ylabel('Outer Planet Migration Speed')
axes[2].set_xlabel(r'$K=\tau_e/\tau_a$')
plt.colorbar(sc, ax=axes[0], label='Period Ratio Amplitude')
plt.colorbar(sc2, ax=axes[1], label='e1 Amplitude')
plt.colorbar(sc3, ax=axes[2], label='e2 Amplitude')
axes[0].set_title('Variable Amplitudes: V$_{max}$ - V$_{mean}$')
'''
