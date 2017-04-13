#This calculates the chi2 values for all 3 different models - PPI, noPPI, L&P (best one)
import numpy as np
import rebound
import pandas as pd
import matplotlib.pyplot as plt

def fit_RV(times,theta,pp_interaction=1):
    m1sini,m2sini,a1,a2,h1,h2,k1,k2,lambda1,lambda2,sini = theta
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
    mJ = 9.543e-4                        #Jupiter mass -> solar mass
    v = np.empty(0)
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.add(m=0.92)                      #add the star
    if pp_interaction == 0:
        sim.testparticle_type = 1
        sim.N_active = sim.N
    sim.add(m=m1sini*mJ/sini,a=a1,l=lambda1,h=h1,k=k1)
    sim.add(m=m2sini*mJ/sini,a=a2,l=lambda2,h=h2,k=k2)
    sim.dt = 2*np.pi* a1**(1.5) / 200.
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
y, yerr2 = data["RV"], data["Unc"]**2

#PPI
samples_PPI = np.load('emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy')[:,1000:,:].reshape((-1, 13))
MAP_PPI = np.percentile(samples_PPI, 50, axis=0)
jitter2, offset = MAP_PPI[-2], MAP_PPI[-1]
model = fit_RV(times, MAP_PPI[:11],1) + offset
#Chi2_PPI = 0.5*np.sum((y - model)**2 /(yerr2 + jitter2) + np.log(yerr2 + jitter2)) / 13
Chi2_PPI = np.sum((y - model)**2 /yerr2 ) / (len(model) - 13)

#noPPI
samples_noPPI = np.load('emcee_chains/best_runs/hksemi-active_400walk_5000it_chkpt1.npy')[:,1000:,:].reshape((-1, 13))
MAP_noPPI = np.percentile(samples_noPPI, 50, axis=0)
jitter2, offset = MAP_noPPI[-2], MAP_noPPI[-1]
model = fit_RV(times, MAP_noPPI[:11],0) + offset
#Chi2_noPPI = 0.5*np.sum((y - model)**2 /(yerr2 + jitter2) + np.log(yerr2 + jitter2)) / 13
Chi2_noPPI = np.sum((y - model)**2 /(yerr2 )) / (len(model) - 13)

#L&P
dtoyr2pi = 2*np.pi/365.
filename='Nbody/saved_output/round23_GOOD_ONES_FINAL/taueinner_migrate3.0e+03_Kin6.7e+00_Kout2.7e+02_sd931'
samples_LP = np.load('%s_chain.npy'%filename)[:,500:,:].reshape((-1, 6))
MAP_LP = np.percentile(samples_LP, 50, axis=0)
x_s, x_t, y_s, y_t, phi, jitter2 = MAP_LP
times = (times*dtoyr2pi + x_t)/x_s
model = y_s*get_simRV(filename,times,phi) + y_t
#Chi2_LP = 0.5*np.sum((y - model)**2/(yerr2 + jitter2) + np.log(yerr2 + jitter2)) / 8
Chi2_LP = 0.5*np.sum((y - model)**2/(yerr2)) / (len(model) - 8)

#plt.errorbar(times, (y - model)/data["Unc"], yerr=data["Unc"], fmt='o', color='green')
#plt.plot([0,60],[0,0],'k--')
#plt.show()

print "Chi^2_reduced:"
print "PPI=%f, noPPI=%f, L&P=%f"%(Chi2_PPI, Chi2_noPPI, Chi2_LP)


