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
    sim.dt = 2*np.pi* a1**(1.5) / 25.
    sim.move_to_com()
    for t in times*dtoyr2pi:
        sim.integrate(t)
        v = np.append(v,-AUyr2ms*sim.particles[0].vy*sini)
    return v

#simulated system
name = sys.argv[1]
samples = np.load(name+'.npy')[:,100:,:].reshape((-1, 3))
sim_MAP = np.percentile(samples, 50, axis=0)
sim_time, sim_RVx, sim_RVy, dE, N = np.loadtxt(open(name+'_RV.txt', 'r'), delimiter=',', unpack=True)

#get Simulated MAP (migration run)
sim_time = sim_MAP[0]*(sim_time + sim_MAP[1])
sim_RV = sim_RVx*np.sin(sim_MAP[2]) + sim_RVy*np.cos(sim_MAP[2])

#MCMC MAP
MAP = np.percentile(np.load('../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy')[:,1000:,:].reshape((-1, 13)), 50, axis=0)[:-2]
MAP_RV = get_RV(sim_time, MAP)

#compare
plt.plot(sim_time, MAP_RV)
plt.plot(sim_time, sim_RV)
plt.savefig(name+'.png')

plot_corner = 1
if plot_corner == 1:
    fig = corner.corner(samples, labels=["x_s", "x_t", "phi"])
    fig.savefig(name+"_corner.png")
