#This is a compliment to fit_sims_to_RV.py, where the emcee run has finished, and now I'm creating a corner plot and seeing how well the RV curve fits the MAP value.

import sys
import rebound
import corner
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import glob
import os
from matplotlib import gridspec

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

fontsize=13

#main code
dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
files = glob.glob(dir+'*.npy')
data = pd.read_csv('../RV.txt', delimiter=' ')
time_RVdays, data_RV, err_RV = data['BJD'], data['RV'], data['Unc']
time_RV = (time_RVdays - time_RVdays[0])*dtoyr2pi
print "analyzing %d files"%(len(files))
for i,f in enumerate(files):
    name = f.split('.npy')[0]
    sim_samples = np.load(name+'.npy')[:,500:,:].reshape((-1, n_params))
    sim_MAP = np.percentile(sim_samples, 50, axis=0)
    x_s, x_t, y_s, y_t, phi, jitter2 = sim_MAP

    #plot sim, data, full curve - plot 1
    fig = plt.figure(figsize=(11,6))
    gs = gridspec.GridSpec(2, 1, height_ratios=[2, 1, 1])
    ax0 = plt.subplot(gs[0])
    ax1 = plt.subplot(gs[1], sharex=ax0)
    plt.subplots_adjust(hspace = 0.25)
    
    times = (time_RV + x_t)/x_s
    simRV = y_s*get_simRV(name,times,phi) + y_t
    times = np.linspace(0,80/x_s,300)                                                   #plot full curve
    simRVfull = y_s*get_simRV(name,times,phi) + y_t
    timesfull = (times*x_s - x_t)/dtoyr2pi + data['BJD'][0]
    lnL = 0.5*np.sum( (simRV - data_RV)**2/(err_RV**2 + jitter2) + np.log(err_RV**2 + jitter2) )
    if lnL < 500 and jitter2 < 20:
        for theta_sample in sim_samples[np.random.randint(len(sim_samples), size=100)]:
            x_s, x_t, y_s, y_t, phi, jitter2 = theta_sample
            times_sample = np.linspace(0,80/x_s,300)
            simRVfull_sample = y_s*get_simRV(name,times_sample,phi) + y_t
            timesfull_sample = (times*x_s - x_t)/dtoyr2pi + data['BJD'][0]
            index = (timesfull_sample<time_RVdays.iloc[-1]) & (timesfull_sample>time_RVdays.iloc[0])
            ax0.plot(timesfull_sample[index],simRVfull_sample[index],color="k", alpha=0.2)
        ax0.plot(timesfull_sample[index],simRVfull_sample[index],color="k", label='posterior samples')

    index = (timesfull<time_RVdays.iloc[-1]+1) & (timesfull>time_RVdays.iloc[0]-1)
    ax0.plot(timesfull[index],simRVfull[index], color='green',linewidth=2, label='MAP curve')
    ax0.errorbar(time_RVdays,data_RV, yerr=err_RV, fmt='.', color='blue', label='data')     #data
    ax0.plot(time_RVdays, simRV, '.', color='red', label='MAP points')                      #sims
    ax0.legend(loc='upper right',fontsize=6,numpoints=1)
    ax0.set_ylabel('RV (m/s)',fontsize=fontsize)
    ax0.set_xlabel('BJD - 2450000',fontsize=fontsize)
    ax0.set_xlim([2000,6100])
    #axes.set_title('lnL = %f'%lnL)
    
    ax1.errorbar(time_RVdays, simRV - data_RV, yerr=err_RV, fmt='.', color='green')
    ax1.plot([time_RVdays.iloc[0],time_RVdays.iloc[-1]+1], [0,0],'k--', lw=2)
    ax1.set_ylabel('MAP Residuals (m/s)',fontsize=fontsize)
    ax1.set_xlabel('BJD - 2450000',fontsize=fontsize)
    if lnL < 500 and jitter2 < 20:
        plt.savefig("%s.pdf"%name)
        fig = corner.corner(sim_samples, labels=["x_s", "x_t", "y_s", "y_t", "phi", "jitter2"])
        fig.savefig("%s_corner.png"%name)
        plt.close(fig)
        os.system("python orbits.py %s.txt"%name)
        dir = name.split('/')
        #os.system("mv %s* %s/%s/good_ones/."%(name,dir[0],dir[1]))
    else:
        plt.savefig("%s.png"%name)
    plt.close()

    if plot_corner == 1:
        fig = corner.corner(sim_samples, labels=["x_s", "x_t", "y_s", "y_t", "phi", "jitter2"])
        fig.savefig(name+"_corner.png")
        plt.close(fig)

    print i
