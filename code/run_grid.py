import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

from scipy.optimize import minimize
import time

import emcee
from gatspy.periodic import LombScargleFast
import corner

import prob
import orbit
import gaia_tools as gaia
import const as c


def load_data(file_name):

    data = np.genfromtxt(file_name)

    names = ["time","ra","dec","ra_err","dec_err","rv","rv_err","plx","plx_err"]
    obs_pos = np.recarray(len(data[:,0]), names=names,
                          formats=['float64,float64,float64,float64,float64,float64,float64,float64,float64'])

    obs_pos['time'] = data[:,0]  # in seconds
    obs_pos['ra'] = data[:,1]/3600.0/1.0e3  # in degrees
    obs_pos['dec'] = data[:,2]/3600.0/1.0e3  # in degrees
    obs_pos['rv'] = (2./3.) * data[:,3]/1.0e3  # in km/s

    obs_pos['ra_err'] = data[:,4]/1.0e6 # in asec
    obs_pos['dec_err'] = data[:,4]/1.0e6 # in asec
    obs_pos['rv_err'] = 1.0

    return obs_pos




def adjust_distance(obs_pos, distance_new, distance_current=2003.0, mag_current=12.5):
    """Adjust object according to new distance

    """

    V_IC = -0.5417

    mag_new = gaia.calc_new_magnitude(mag_current, distance_current, distance_new)

    plx_err = gaia.get_plx_err(G=mag_new, V_IC=V_IC) * 1.0e3  # plx error in mas

    obs_pos['plx'] = 1.0e3 / distance_new
    obs_pos['plx_err'] = plx_err


    # Adjust position coordinates
    obs_pos['ra'] *= distance_current / distance_new
    obs_pos['dec'] *= distance_current / distance_new

    obs_pos['ra_err'] = gaia.get_single_obs_err(G=mag_new, V_IC=V_IC)
    obs_pos['dec_err'] = gaia.get_single_obs_err(G=mag_new, V_IC=V_IC)

    return obs_pos

def add_uncertainties(obs_pos):
    obs_pos['ra'] += np.random.normal(scale=obs_pos['ra_err']/3600.0, size=len(obs_pos))
    obs_pos['dec'] += np.random.normal(scale=obs_pos['dec_err']/3600.0, size=len(obs_pos))
    obs_pos['rv'] += np.random.normal(scale=obs_pos['rv_err'], size=len(obs_pos))

    return obs_pos


def get_truths(distance_new):

    sys_ra = 1.0e-10
    sys_dec = 1.0e-10
    Omega = 2.1
    omega = 2.5
    I = 0.3
    # tau = 0.33 + P/2.1
    e = 0.7
    P = 2.0*3600.0*24.0*365.25
    tau = P/500.0
    tau = 100.0
    gamma = 1.0e-4
    M1 = 10.0*c.Msun
    M2 = 20.0*c.Msun

    truths = sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance_new

    return truths


def initialize_sampler(truths, obs_pos, nwalkers=100):

    ndim = 12

    # Starting array
    p0 = np.zeros((nwalkers, ndim))

    # Random deviations from best fit
    delta = 1.0 + np.random.normal(0.0, 0.001, size=(nwalkers, ndim))

    # Adjust starting array accordingly
    for i in range(ndim): p0[:,i] = truths[i] * delta[:,i]

    return p0


def run_emcee(p0, obs_pos):

    nwalkers, ndim = p0.shape

    # Define sampler
    sampler = emcee.EnsembleSampler(nwalkers=nwalkers, dim=ndim, lnpostfn=prob.ln_posterior, args=(obs_pos,))

    # Run burn-in
    pos,ln_prob,state = sampler.run_mcmc(p0, N=1000)

    # Run sampler
    sampler.reset()
    pos,ln_prob,state = sampler.run_mcmc(pos, N=5000)

    # Downsample by factor of 10
    chains = sampler.chain[:,::10,:]

    return chains


def run_one_binary(distance, rv_err):

    # Load binary data
    obs_pos = load_data("../data/Eccentric_True_RA_Dec_RV.txt")

    # Adjust binary to new distance
    obs_pos = adjust_distance(obs_pos, distance)

    # Adjust binary's RV error
    obs_pos['rv_err'] = rv_err

    # Add uncertainties to position and RV measurements
    obs_pos = add_uncertainties(obs_pos)

    # Calculate truths for the binary
    truths = get_truths(distance)

    # Initialize the emcee sampler
    p0 = initialize_sampler(truths, obs_pos, nwalkers=100)

    # Run the sampler
    chains = run_emcee(p0, obs_pos)

    # Save the chains
    file_out = "../data/binary_chains/distance" + "%.2f" % distance + "_rverr" + "%.2f" % rv_err + ".npy"
    np.save(file_out, chains)



def run_grid(n_dist=10, n_rv_err=10, dist_range=(50, 10000), rv_err_range=(0.3, 20.0)):

    start_time = time.time()

    dist_set = 10**np.linspace(np.log10(dist_range[0]), np.log10(dist_range[1]), n_dist)
    rv_err_set = 10**np.linspace(np.log10(rv_err_range[0]), np.log10(rv_err_range[1]), n_rv_err)

    for dist in dist_set:
        for rv_err in rv_err_set:

            print("Running:", dist, rv_err, time.time()-start_time)
            run_one_binary(dist, rv_err)


run_grid()
