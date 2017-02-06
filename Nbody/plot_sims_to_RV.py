#This is a compliment to fit_sims_to_RV.py, where the emcee run has finished, and now I'm creating a corner plot and seeing how well the RV curve fits the MAP value.

import sys
import rebound
import corner
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

def get_simRV(filename, time_sim, phi):
    dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    sim = rebound.Simulation.from_file(filename+'.bin')
    tmax = sim.t
    rv = np.empty(0)
    for t in time_sim*dtoyr2pi:
        sim.integrate(t+tmax,1)
        rv = np.append(rv,AUyr2ms*(sim.particles[0].vx*np.sin(phi) + sim.particles[0].vy*np.cos(phi)))
    del sim
    return rv

#args
dir = sys.argv[1]
plot_corner = 1

#main code
files = glob.glob(dir+'*.npy')
data = pd.read_csv('../RV.txt', delimiter=' ')
for f in files:
    name = f.split('.npy')[0]
    sim_samples = np.load(name+'.npy')[:,500:,:].reshape((-1, 5))
    sim_MAP = np.percentile(sim_samples, 50, axis=0)
    x_s, x_t, phi, jitter2, offset = sim_MAP

    #plot data
    t = data['BJD']-data['BJD'][0]
    plt.errorbar(t, data['RV'], yerr=data['Unc'], fmt='o')

    #plot sim
    times = np.linspace(0,(data['BJD'].iloc[-1]-data['BJD'].iloc[0] + 50 - x_t)/x_s,200)
    simRV = get_simRV(name,times,phi) + offset
    t_sim = times*x_s + x_t
    plt.plot(t_sim, simRV, label='sim', color='red')

    plt.legend(loc='upper left', fontsize=8, numpoints=1)
    plt.savefig(name+".png")
    plt.clf()

    if plot_corner == 1:
        fig = corner.corner(sim_samples, labels=["x_s", "x_t", "phi", "jitter2", "offset"])
        fig.savefig(name+"_corner.png")
        plt.close(fig)
