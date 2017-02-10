#This is a compliment to fit_sims_to_RV.py, where the emcee run has finished, and now I'm creating a corner plot and seeing how well the RV curve fits the MAP value.

import sys
import rebound
import corner
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob

def get_simRV(filename, time_sim, phi):
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    sim = rebound.Simulation.from_file(filename+'.bin')
    tmax = sim.t
    rv = np.empty(0)
    for t in time_sim:
        sim.integrate(t+tmax,1)
        rv = np.append(rv,AUyr2ms*(sim.particles[0].vx*np.sin(phi) + sim.particles[0].vy*np.cos(phi)))
    del sim
    return rv

#args
dir = sys.argv[1]
plot_corner = 0
n_params = 6    #x_s, x_t, y_s, y_t, phi, jitter2

#main code
dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
files = glob.glob(dir+'*.npy')
data = pd.read_csv('../RV.txt', delimiter=' ')
time_RV, data_RV, err_RV = (data['BJD']-data['BJD'][0])*dtoyr2pi, data['RV'], data['Unc']
print "analyzing %d files"%(len(files))
for i,f in enumerate(files):
    fig = plt.figure(figsize=(11,6))
    name = f.split('.npy')[0]
    sim_samples = np.load(name+'.npy')[:,500:,:].reshape((-1, n_params))
    sim_MAP = np.percentile(sim_samples, 50, axis=0)
    x_s, x_t, y_s, y_t, phi, jitter2 = sim_MAP

    #plot sim, data, full curve - plot 1
    axes = fig.add_subplot(2, 1, 1)
    times = (time_RV + x_t)/x_s
    simRV = y_s*get_simRV(name,times,phi) + y_t
    axes.errorbar(time_RV, data_RV, yerr=err_RV, fmt='.', label='data points')          #data
    axes.plot(time_RV, simRV, '.', color='red', label='sim points')                     #sims
    times = np.linspace(0,80/x_s,300)                                                   #plot full curve
    simRVfull = y_s*get_simRV(name,times,phi) + y_t
    timesfull = times*x_s - x_t
    axes.plot(timesfull[timesfull<(time_RV.iloc[-1]+1)], simRVfull[timesfull<(time_RV.iloc[-1]+1)], color='red', lw=0.5, label='full sim curve')
    axes.legend(loc='lower left', fontsize=8, numpoints=1)
    lnL = 0.5*np.sum( (simRV - data_RV)**2/(err_RV**2 + jitter2) + np.log(err_RV**2 + jitter2) )
    axes.set_ylabel('RV (m/s)')
    axes.set_title('lnL = %f'%lnL)
    
    axes = fig.add_subplot(2, 1, 2, sharex=axes)
    axes.errorbar(time_RV, simRV - data_RV, yerr=err_RV, fmt='.')
    axes.plot([time_RV.iloc[0],time_RV.iloc[-1]+1], [0,0], 'k--')
    axes.set_xlabel('time')
    axes.set_ylabel('residuals (m/s)')
    if lnL < 500 and jitter2 < 20:
        dir = name.split('/')
        plt.savefig("%s/%s/good_ones/%s.png"%(dir[0],dir[1],dir[2]))
        fig = corner.corner(sim_samples, labels=["x_s", "x_t", "y_s", "y_t", "phi", "jitter2"])
        fig.savefig("%s/%s/good_ones/%s_corner.png"%(dir[0],dir[1],dir[2]))
        plt.close(fig)
    else:
        plt.savefig("%s.png"%name)
    plt.close()

    if plot_corner == 1:
        fig = corner.corner(sim_samples, labels=["x_s", "x_t", "y_s", "y_t", "phi", "jitter2"])
        fig.savefig(name+"_corner.png")
        plt.close(fig)

    print i
