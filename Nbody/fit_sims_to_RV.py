#This script:
#a) Draws masses from the MCMC posterior chain, runs migration simulations (K1, K2, mig_rate = free parameters). Run the C version of this code (faster, easier). 
#b) After migration ends, get some RV outputs and run an MCMC chain comparing it to the original data (x-stretch, x-translate, viewing angle = free parameters).
#c) See if you can find optimal parameters this way.

import multiprocessing as mp
import os
import sys
import random
import numpy as np
import pandas as pd
import emcee
import rebound
from progress.bar import Bar

#############General Functions######################
def make_runs(N_runs):
    #draw masses from the posterior
    m1 = []
    m2 = []
    sini = []
    burnin = 1000
    ndim = 13
    filename = '../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy'
    samples = np.load(filename)[:, burnin:, :].reshape((-1, ndim))
    for theta in samples[np.random.randint(len(samples), size=N_runs)]:
        m1.append(theta[0])
        m2.append(theta[1])
        sini.append(theta[10])
    #make N_runs for simulation
    random.seed()
    runs = []
    mig_rate = random.sample(np.round(np.logspace(2,4,N_runs)), N_runs)
    K1 = random.sample(np.logspace(-1,2,N_runs), N_runs)
    K2 = np.logspace(-1,2,N_runs)
    path = 'output/'
    for i in xrange(0,N_runs):
        seed = int(1000*random.random())
        name = path+'taueinner_migrate%.1e_Kin%.1e_Kout%.1e_sd%d'%(mig_rate[i],K1[i],K2[i],seed)
        runs.append((m1[i],m2[i],sini[i],mig_rate[i],K1[i],K2[i],seed,name))
    return runs

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

#############emcee stuff############################
def lnlike(theta, filename, time_RV, data_RV, err2_RV):
    x_s, x_t, phi, jitter2, offset = theta
    time_sim = (time_RV - x_t)/x_s
    sim_RV = get_simRV(filename, time_sim, phi) + offset
    return -0.5*np.sum( (sim_RV - data_RV)**2/(err2_RV + jitter2) + np.log(err2_RV + jitter2) )

def lnprior(theta):
    x_s, x_t, phi, jitter2, offset = theta        #x-stretch, x-translate, phi (viewing angle), jitter2, RV offset
    if 0.5<x_s<2.5 and -600<x_t<0 and 0<=phi<2*np.pi and 0.<jitter2<500. and -40<offset<40:
        return 0
    return -np.inf

def lnprob(theta, filename, time_RV, data_RV, err2_RV):
    lnp = lnprior(theta)
    if not np.isfinite(lnp):
        return -np.inf
    lnL = lnlike(theta, filename, time_RV, data_RV, err2_RV)
#    f = open('output/emcee_out.txt','a')
#    f.write('x_s=%f \t x_t=%f \t phi=%f \t j2=%f \t off=%f \t lnp=%f, lnL=%f\n'%(theta[0],theta[1],theta[2],theta[3],theta[4],lnp, lnL))
#    f.close()
    return lnp + lnL

def run_emcee(filename, time_RV, data_RV, err2_RV):
    theta_ini = [1.5,0,np.pi,10,4]
    ndim, nwalkers, n_it, bar_checkpoints = len(theta_ini), 26, 2000, 100
    pos = [theta_ini + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(filename, time_RV, data_RV, err2_RV));
    bar = Bar('Processing', max=bar_checkpoints)
    for i in range(bar_checkpoints):
        pos, _, _ = sampler.run_mcmc(pos, n_it/bar_checkpoints);
        bar.next()
    bar.finish()
    np.save(filename+".npy",sampler.chain)

####################################################
#############Main Code##############################
samples = np.load('../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy')[:,1000:,:].reshape((-1, 13))
MAPP = np.percentile(samples, 50, axis=0)[:-2]

#each pool worker executes this
def execute(pars):
    os.system('./rebound %f %f %f %f %f %f %d %s'%pars)
    name = pars[-1].split('.txt')[0]
    try:
        print "\nPerforming MCMC fit."
        data = pd.read_csv('../RV.txt', delimiter=' ')
        time_RV, data_RV, err2_RV = data['BJD'] - data['BJD'][0], data['RV'], data['Unc']**2
        run_emcee(name, time_RV, data_RV, err2_RV)
    except:
        f = open("output/bad_sims.txt","a")
        f.write("Error simulating %s.txt. Skipped emcee.\n"%name)
        f.close()
        print "\nError simulating %s.txt. Skipping emcee.\n"%name

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    os.system('make')
    N_runs = 20
    pool = mp.Pool(processes=np.min([N_runs, 2]))
    runs = make_runs(N_runs)
    #runs = [(0.90721388757667032, 0.8489328864365624, 0.95085548551813603, 2000.0, 1.0, 0.1, 646, 'output/taueinner_migrate1.0e+03_Kin1.0_Kout1.0_sd646')]
    #runs = [(0.99672170557149731, 0.87038713759372832, 0.82730589001482202, 1000.0, 5.0, 5.0, 757, 'output/taueinner_migrate1.0e+03_Kin5.0e+00_Kout5.0e+00_sd757')]
    pool.map(execute, runs)
    pool.close()
    pool.join()
