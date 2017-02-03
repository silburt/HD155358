#This is a compliment to fit_sims_to_RV.py, where the emcee run has finished, and now I'm creating a corner plot and seeing how well the RV curve fits the MAP value.

import sys
import rebound
import corner
import numpy as np
import matplotlib.pyplot as plt

def get_RV(times,theta):
    m1sini,m2sini,a1,a2,h1,h2,k1,k2,lambda1,lambda2,sini = theta
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
    mJ = 9.543e-4                        #Jupiter mass -> solar mass
    v = np.empty(0)
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.add(m=0.92)
    sim.add(m=m1sini*mJ/sini,a=a1,l=lambda1,h=h1,k=k1)
    sim.add(m=m2sini*mJ/sini,a=a2,l=lambda2,h=h2,k=k2)
    sim.dt = 2*np.pi* a1**(1.5) / 31.
    sim.move_to_com()
    for t in times*dtoyr2pi:
        sim.integrate(t, 1)
        v = np.append(v,-AUyr2ms*sim.particles[0].vy*sini)
    return v

#simulated system
name = sys.argv[1]
#samples = np.load(name+'.npy')[:,100:,:].reshape((-1, 3))
#sim_MAP = np.percentile(samples, 50, axis=0)
#x_s, x_t, phi = sim_MAP
time, sim_RVx, sim_RVy, dE, N = np.loadtxt(open(name+'_RV.txt', 'r'), delimiter=',', unpack=True)

x_s, x_t, phi = 1, 0, 3.03

#MCMC MAP
MAP = np.percentile(np.load('../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy')[:,1000:,:].reshape((-1, 13)), 50, axis=0)[:-2]

vals = [(x, y) for x in np.linspace(0,2*np.pi,10) for y in np.linspace(1,3,10)]

lnL = []
for phi,xs in vals[0:1]:
    sim_time = xs*(time)
    MAP_RV = get_RV(sim_time, MAP)
    sim_RV = sim_RVx*np.sin(phi) + sim_RVy*np.cos(phi)
    lnL.append(-0.5*np.sum( (MAP_RV - sim_RV)**2 ))
    print xs,phi



'''
#compare
for x_s in np.linspace(1.45,1.55,30):
    sim_time = x_s*(time + x_t)
    MAP_RV = get_RV(sim_time, MAP)
    sim_RV = sim_RVx*np.sin(phi) + sim_RVy*np.cos(phi)
    plt.plot(sim_time/sim_time[-1], MAP_RV)
    plt.plot(sim_time/sim_time[-1], sim_RV)
    plt.title('lnL=%f'%(-0.5*np.sum( (MAP_RV - sim_RV)**2 )))
    plt.savefig('output/phi%.3f.png'%x_s)
    plt.clf()
    print "hi"
'''

plot_corner = 0
if plot_corner == 1:
    fig = corner.corner(samples, labels=["x_s", "x_t", "phi"])
    fig.savefig(name+"_corner.png")