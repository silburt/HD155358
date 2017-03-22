#This calculates the chi2 values for all 3 different models - PPI, noPPI, L&P (best one)
import numpy as np
import rebound
import pandas as pd

def fit_RV(times,theta,pp_interaction=1):
    m1sini,m2sini,a1,a2,h1,h2,k1,k2,lambda1,lambda2,sini = theta
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
    mJ = 9.543e-4                        #Jupiter mass -> solar mass
    v = np.empty(0)
    sim = rebound.Simulation()
    sim.integrator = "ias15"
    sim.add(m=0.92)                      #add the star
    if pp_interaction == 0:
        sim.testparticle_type = 1
        sim.N_active = sim.N
    sim.add(m=m1sini*mJ/sini,a=a1,l=lambda1,h=h1,k=k1)
    sim.add(m=m2sini*mJ/sini,a=a2,l=lambda2,h=h2,k=k2)
    sim.dt = 2*np.pi* a1**(1.5) / 25.    #dt = 25 steps per orb per of inner planet
    sim.move_to_com()
    for t in times*dtoyr2pi:
        sim.integrate(t)
        v = np.append(v,-AUyr2ms*sim.particles[0].vy*sini)
    return v

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

#data
data = pd.read_csv("RV.txt", delimiter=' ')
times = data["BJD"] - data["BJD"].iloc[0]

#PPI
samples_PPI = np.load('emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy')[:,1000:,:].reshape((-1, 13))
MAP_PPI = np.percentile(samples_PPI, 50, axis=0)
Chi2_PPI = np.sum((data["RV"] - (fit_RV(times,MAP_PPI[:11],1) + MAP_PPI[12]))**2 / data["Unc"]**2) / 13

#noPPI
samples_noPPI = np.load('emcee_chains/best_runs/hksemi-active_400walk_5000it_chkpt1.npy')[:,1000:,:].reshape((-1, 13))
MAP_noPPI = np.percentile(samples_noPPI, 50, axis=0)
Chi2_noPPI = np.sum((data["RV"] - (fit_RV(times,MAP_noPPI[:11],0) + MAP_noPPI[12]))**2 / data["Unc"]**2) /13

#L&P
dtoyr2pi = 2*np.pi/365.
filename='Nbody/saved_output/round22_GOOD_ONES_aggregated/taueinner_migrate3.0e+03_Kin6.7e+00_Kout2.7e+02_sd931'
samples_LP = np.load('%s_chain.npy'%filename)[:,500:,:].reshape((-1, 6))
MAP_LP = np.percentile(samples_LP, 50, axis=0)
x_s, x_t, y_s, y_t, phi, jitter2 = MAP_LP
times = (times*dtoyr2pi + x_t)/x_s
Chi2_LP = np.sum((data["RV"] - (y_s*get_simRV(filename,times,phi) + y_t))**2 / data["Unc"]**2) / 6

print Chi2_PPI, Chi2_noPPI, Chi2_LP


