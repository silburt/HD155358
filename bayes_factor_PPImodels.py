#This computes the Bayes' factor between the inclination run and the no-inclination run for our PPI model

import numpy as np
import pandas as pd
import sys
import matplotlib.pyplot as plt
import rebound

def fit_RV(times,theta):
    if len(theta) == 13:
        m1sini,m2sini,a1,a2,h1,h2,k1,k2,lambda1,lambda2,sini = theta[:11]
        ix2, iy2 = 0, 0
    else:
        m1sini,m2sini,a1,a2,h1,h2,k1,k2,lambda1,lambda2,ix2,iy2,sini = theta[:13]
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    dtoyr2pi = 2*np.pi/365.              #days -> yr/2pi
    mJ, m0 = 9.543e-4, 0.92              #Jupiter mass -> solar mass, star mass
    v = np.empty(0)
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.add(m=m0)                        #add the star
    sim.add(m=m1sini*mJ/sini,a=a1,l=lambda1,h=h1,k=k1)
    sim.add(m=m2sini*mJ/sini,a=a2,l=lambda2,h=h2,k=k2,ix=ix2,iy=iy2)
    sim.dt = 2*np.pi*a1**(1.5)/m0/200.    #dt = 200 steps per orb per of inner planet
    sim.move_to_com()
    for t in times*dtoyr2pi:
        sim.integrate(t)
        v = np.append(v,-AUyr2ms*sim.particles[0].vy*sini)
    return v

def get_MLE(samples, data, n_params):
    theta = np.percentile(samples, 50, axis=0)
    jitter2, offset = theta[-2], theta[-1]
    times = data["BJD"] - data["BJD"].iloc[0]               #BJD (days)
    y, yerr2 = data["RV"].values, (data["Unc"].values)**2
    model = fit_RV(times, theta) + offset
    return np.exp(-0.5*np.sum( (y - model)**2/(yerr2 + jitter2) + np.log(yerr2 + jitter2) ))

data = pd.read_csv("RV.txt", delimiter=' ')

#PPI model - no inclination
n_params = 13
samples = np.load('emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy')[:,1000:,:].reshape((-1, n_params))
MLE = get_MLE(samples,data,n_params)
std = np.sqrt(np.linalg.det(np.cov(samples, rowvar=0))) * (2*np.pi)**(n_params/2.)
L1 = MLE*std

#PPI model - inclination
n_params = 15
samples = np.load('emcee_chains/best_runs/hk-inc_300walk_6000it.npy')[:,2000:,:].reshape((-1, n_params))
MLE = get_MLE(samples,data,n_params)
std = np.sqrt(np.linalg.det(np.cov(samples, rowvar=0))) * (2*np.pi)**(n_params/2.)
L2 = MLE*std*(1./16)

print L1/L2

