#This script:
#a) Draws masses from the MCMC posterior chain, runs migration simulations (K1, K2, mig_rate = free parameters). Run the C version of this code (faster, easier). 
#b) After migration ends, get some RV outputs and run an MCMC chain comparing it to the MAP_RV value (x-stretch, x-translate, viewing angle = free parameters).
#c) See if you can find optimal parameters this way.

import multiprocessing as mp
import os
import sys
import random
import numpy as np
import emcee
import rebound
from progress.bar import Bar
import peakutils

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

def sim_MAP(times,theta):
    m1sini,m2sini,a1,a2,h1,h2,k1,k2,lambda1,lambda2,sini = theta
    AUyr2ms = 29682.77                   #AU/(yr/2pi) -> m/s
    mJ = 9.543e-4                        #Jupiter mass -> solar mass
    v = np.empty(0)
    sim = rebound.Simulation()
    sim.integrator = "whfast"
    sim.add(m=0.92)
    sim.add(m=m1sini*mJ/sini,a=a1,l=lambda1,h=h1,k=k1)
    sim.add(m=m2sini*mJ/sini,a=a2,l=lambda2,h=h2,k=k2)
    sim.dt = 2*np.pi* a1**(1.5) / 31.
    sim.move_to_com()
    for t in times:
        sim.integrate(t,1)
        v = np.append(v,-AUyr2ms*sim.particles[0].vy*sini)
    return v

#Initial Value stuff (need to get close to the real solution before emcee will work)
#############Initial Value##########################
def get_theta_ini(time, sim_RVx, sim_RVy, MAP):
    print "Getting initial conditions..."
    MAP_RV = sim_MAP(time, MAP)
    
    #stretch factor
    i_sim = peakutils.indexes(sim_RVy, thres=0.9) #locate successive peaks to get period
    i_MAP = peakutils.indexes(-1*MAP_RV, thres=0.9)
    P_sim = time[i_sim[1]] - time[i_sim[0]]
    P_MAP = time[i_MAP[1]] - time[i_MAP[0]]
    xs = P_MAP/P_sim

    #phi
    MSE_phi = []   #mean squred error of max peak values + min peak values
    phi_angles = np.linspace(0,2*np.pi,50)
    for p in phi_angles:
        sim_RV = sim_RVx*np.sin(p) + sim_RVy*np.cos(p)
        MSE_phi.append( (np.max(sim_RV)-np.max(MAP_RV))**2 + (np.min(sim_RV)-np.min(MAP_RV))**2 )
    MSE_phi = np.asarray(MSE_phi)
    phi = phi_angles[MSE_phi == np.min(MSE_phi)][0]

    #translation factor
    xt_vals = np.linspace(-P_MAP/2,P_MAP/2,20, endpoint=False)
    MSE_xt = []
    for tl in xt_vals:
        sim_time = xs*time + tl
        MAP_RV = sim_MAP(sim_time, MAP)
        sim_RV = sim_RVx*np.sin(phi) + sim_RVy*np.cos(phi)
        MSE_xt.append( np.sum((MAP_RV - sim_RV)**2) )
    MSE_xt = np.asarray(MSE_xt)
    xt = xt_vals[MSE_xt == np.min(MSE_xt)][0]
    return (xs, xt, phi)

#############emcee stuff############################
def lnlike(theta, sim_time, sim_RVx, sim_RVy, MAP):
    x_s, x_t, phi = theta
    sim_time = x_s*sim_time + x_t
    sim_RV = sim_RVx*np.sin(phi) + sim_RVy*np.cos(phi)
    MAP_RV = sim_MAP(sim_time, MAP)
    #MAP_err2 = (len(sim_time))**2           #each MAP data point has some uncertainty to it.
    #return -0.5*np.sum( (sim_RV - MAP_RV)**2/MAP_err2 + np.log(MAP_err2) )
    return -0.5*np.sum( (sim_RV - MAP_RV)**2 )

def lnprior(theta):
    x_s, x_t, phi = theta        #x-stretch, x-translate, sinphi (viewing angle)
    if 0.5<x_s<2.5 and -5<x_t<5 and 0<=phi<2*np.pi:
        return 0
    return -np.inf

def lnprob(theta, sim_time, sim_RVx, sim_RVy, MAP):
    lnp = lnprior(theta)
    if not np.isfinite(lnp):
        return -np.inf
    lnL = lnlike(theta, sim_time, sim_RVx, sim_RVy, MAP)
    f = open('output/emcee_out.txt','a')
    f.write('x_s=%f, x_t=%f, phi=%f, lnp=%f, lnL=%f\n'%(theta[0], theta[1], theta[2], lnp, lnL))
    f.close()
    return lnp + lnL

def run_emcee(sim_time, sim_RVx, sim_RVy, MAP, filename):
    theta_ini = get_theta_ini(sim_time, sim_RVx, sim_RVy, MAP) #[1,30,np.pi]   #x_stretch, x_translate, phi (viewing angle)
    ndim, nwalkers, n_it, bar_checkpoints = len(theta_ini), 50, 200, 100
    pos = [theta_ini + 1e-4*np.random.randn(ndim) for i in range(nwalkers)]
    sampler = emcee.EnsembleSampler(nwalkers, ndim, lnprob, args=(sim_time, sim_RVx, sim_RVy, MAP));
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
        sim_time, sim_RVx, sim_RVy, dE, N = np.loadtxt(open(name+'_RV.txt', 'r'), delimiter=',', unpack=True)
        print "\nPerforming MCMC fit."
        run_emcee(sim_time, sim_RVx, sim_RVy, MAPP, name)
    except:
        f = open("output/bad_sims.txt","a")
        f.write("Error simulating %s.txt. Skipped emcee.\n"%name)
        f.close()
        print "\nError simulating %s.txt. Skipping emcee.\n"%name

#Main multiprocess execution - Give sysname and letters of outer planets close to resonance
if __name__== '__main__':
    os.system('make')
    N_runs = 1
    pool = mp.Pool(processes=np.min([N_runs, 1]))
    #runs = make_runs(N_runs)
    runs = [(0.99672170557149731, 0.87038713759372832, 0.82730589001482202, 1000.0, 5.0, 5.0, 757, 'output/taueinner_migrate1.0e+03_Kin5.0e+00_Kout5.0e+00_sd757')]
    pool.map(execute, runs)
    pool.close()
    pool.join()
