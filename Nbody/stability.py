import rebound
import numpy as np
import os
import random
import multiprocessing as mp
import time

def exit_condition(sim):
    o1 = sim.particles[1].calculate_orbit(sim.particles[0])
    o2 = sim.particles[2].calculate_orbit(sim.particles[0])
    hill_1 = o1.a*(0.3333*sim.particles[1].m/sim.particles[0].m)**(1./3.)
    hill_2 = o2.a*(0.3333*sim.particles[2].m/sim.particles[0].m)**(1./3.)
    return np.min([hill_1,hill_2])

def output(sim,E0,f):
    dE = (E0 - sim.calculate_energy())/E0
    o1 = sim.particles[1].calculate_orbit(sim.particles[0])
    o2 = sim.particles[2].calculate_orbit(sim.particles[0])
    f.write("%e,%e,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n"%(sim.t/(2*np.pi),dE,o1.a,o1.e,o1.l,o1.pomega,o1.M,o2.a,o2.e,o2.l,o2.pomega,o2.M))

def simulate(theta,times,name):
    m1sini,m2sini,aa1,aa2,h1,h2,k1,k2,lambda1,lambda2,sini = theta
    mJ = 9.543e-4                       #Jupiter mass -> solar mass
    
    sim = rebound.Simulation()
    sim.integrator = 'whfast'
    sim.dt = 2*np.pi*aa1**(1.5) / 200.  #dt = 100 steps per orb per of inner planet
    sim.add(m=0.92)                     #add the star
    sim.add(m=m1sini*mJ/sini,a=aa1,l=lambda1,h=h1,k=k1)
    sim.add(m=m2sini*mJ/sini,a=aa2,l=lambda2,h=h2,k=k2)
    sim.move_to_com()
    sim.exit_min_distance = exit_condition(sim)
    E0 = sim.calculate_energy()
    with open(name, 'w') as f:
        output(sim,E0,f)
        try:
            for t in times:
                sim.integrate(t)
                output(sim,E0,f)
        except rebound.Encounter as error:
            output(sim,E0,f)
    return sim.t/(2*np.pi)

def make_runs(N_runs, logtmax):
    random.seed()
    runs = []
    filename = '../emcee_chains/best_runs/hk_400walk_5000it_chkpt1.npy'
    #filename = '../emcee_chains/best_runs/hk_250walk_6000it/hk_250walk_6000it_chkpt5.npy'
    #filename = '../emcee_chains/best_runs/hksemi-active_400walk_5000it_chkpt1.npy'         #semi-active!
    samples = np.load(filename)[:, 1000:, :].reshape((-1, 13))
    seeds = np.random.randint(len(samples), size=N_runs)
    for i in range(N_runs):
        theta = samples[seeds[i]]
        name = "output/tmax1e%d_sd%d.csv"%(logtmax,seeds[i])
        runs.append((theta[:11],logtmax,name))
    return runs

def execute(pars):
    theta, logtmax, name = pars
    Npoints = 2000
    times = 2*np.pi*np.logspace(1,logtmax,Npoints)
    t = simulate(theta,times,name)
    print "%e simulated years for %s"%(t,name)

#Main multiprocess execution
if __name__== '__main__':
    #params
    N_runs = 500
    logtmax = 9

    #run
    runs = make_runs(N_runs,logtmax)
    pool = mp.Pool(processes=np.min([N_runs, 20]))
    pool.map(execute, runs)
    pool.close()
    pool.join()
