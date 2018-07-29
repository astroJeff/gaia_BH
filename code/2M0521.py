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

from run_set_pygaia import create_data_array

import astropy.coordinates as coord
import astropy.units as u



# Mass error on the giant star
M2 = 3.0*c.Msun
M2_err = 1.0*c.Msun

c.time_max = 5.0
c.N_samples = int(213. * 60./22.)



def create_binary(M_BH, rv_err=1.0):

    truths = get_truths_2M0521(M_BH)
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = truths
    parallax = 1.0/distance



    #### Use astropy to get 3D coordinates ####
    coor = coord.SkyCoord(ra=sys_ra*u.degree, dec=sys_dec*u.degree, distance=distance*u.pc)
    gc1 = coor.transform_to(coord.Galactocentric)




    # Create array for 2M0521
    dtype = [("M1","<f8"), ("M2","<f8"), ("ecc","<f8"), ("P_orb","<f8"),
                 ("Lum2", "<f8"), ("Temp2", "<f8"), ("Rad2", "<f8"),
                 ("Xgx", "<f8"), ("Ygx", "<f8"), ("Zgx", "<f8"), ("dist","<f8"),
                 ("inc","<f8"), ("Omega","<f8"), ("omega","<f8"),
                 ("proj_sep","<f8"), ("G","<f8"), ("V_IC","<f8")]

    sys = np.zeros(1, dtype=dtype)


    sys['M1'] = M1  # in Msun
    sys['M2'] = M2  # in Msun
    sys['ecc'] = e
    sys['P_orb'] = P # in years
    sys['dist'] = 1.0/parallax  # in pc
    sys['G'] = 12.27
    sys['V_IC'] = 1.86  # Really this is BP-RP

    sys['Xgx'] = gc1.x/u.pc
    sys['Ygx'] = gc1.y/u.pc
    sys['Zgx'] = gc1.z/u.pc

    # Randomly chosen
    sys['inc'] = I
    sys['Omega'] = Omega
    sys['omega'] = omega


    # Create data
    obs_pos = create_data_array(sys, truths, rv_err=1.0)

    return obs_pos





def add_uncertainties(obs_pos):
    # obs_pos['ra'] += np.random.normal(scale=obs_pos['ra_err']/3600.0, size=len(obs_pos))
    # obs_pos['dec'] += np.random.normal(scale=obs_pos['dec_err']/3600.0, size=len(obs_pos))
    # obs_pos['rv'] += np.random.normal(scale=obs_pos['rv_err'], size=len(obs_pos))

    return obs_pos


def get_truths_2M0521(M_BH):

    # Astrometry from Gaia DR2
    sys_ra = 80.4858333
    sys_dec = 43.9894444
    parallax = 0.27
    distance = 1.0e3/parallax

    # Random values
    Omega = 2.1
    omega = 2.5
    I = 0.3
    # tau = 0.33 + P/2.1

    # Fits from Thompson et al
    e = 0.01
    P = 83.2 * 3600.0*24.0  # 83.2 day orbital period
    tau = 100.0
    gamma = 2.56
    M1 = 3.0*c.Msun
    M2 = M_BH*c.Msun

    truths = sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance

    return truths


def initialize_sampler(truths, obs_pos, nwalkers=100):

    ndim = 12

    # Starting array
    p0 = np.zeros((nwalkers, ndim))

    # Random deviations from best fit
    delta = 1.0 + np.random.normal(0.0, 1.0e-5, size=(nwalkers, ndim))
    delta[:,0] = 1.0 + np.random.normal(0.0, 1.0e-10, size=(nwalkers))
    delta[:,1] = 1.0 + np.random.normal(0.0, 1.0e-10, size=(nwalkers))

    # Adjust starting array accordingly
    for i in range(ndim): p0[:,i] = truths[i] * delta[:,i]

    return p0


def run_emcee(p0, obs_pos, M2, M2_err):

    nwalkers, ndim = p0.shape

    # Define sampler
    sampler = emcee.EnsembleSampler(nwalkers=nwalkers, dim=ndim, lnpostfn=prob.ln_posterior, args=(obs_pos, M2, M2_err))

    # Run burn-in
    pos,ln_prob,state = sampler.run_mcmc(p0, N=1000)

    # Run sampler
    sampler.reset()
    pos,ln_prob,state = sampler.run_mcmc(pos, N=2000)

    # Downsample by factor of 10
    chains = sampler.chain[:,::10,:]

    return chains


def run_one_binary(M_BH, rv_err=1.0):

    # Create new data array - set RV error to 1 km/s
    obs_pos = create_binary(M_BH, rv_err=1.0)

    # Add uncertainties to position and RV measurements
    obs_pos = add_uncertainties(obs_pos)

    np.save("../data/2M0521_obs_pos.npy", obs_pos) 

    # Calculate truths for the binary
    truths = get_truths_2M0521(M_BH)

    # Initialize the emcee sampler
    print("Initializing the sampler...")
    p0 = initialize_sampler(truths, obs_pos, nwalkers=100)
    np.save("../data/2M0521_initialized_sampler_pos.npy", p0)
    print("...initialization complete.")

    # Test initialization
    prob_post = prob.ln_posterior(p0[0], obs_pos, M_BH, 0.5)
    print("Posterior probability of initial position. M_BH = ", M_BH, ": ", prob_post)  

    # Run the sampler
    chains = run_emcee(p0, obs_pos, M2, M2_err)

    # Save the chains
    file_out = "../data/binary_chains/2M0521_" + str(M_BH) + "M.npy"
    np.save(file_out, chains)



def run_grid():

    start_time = time.time()

    M_set = [1.4, 2.0, 4.0, 6.0, 8.0]

    for mass in M_set:

        print("Running:", mass, time.time()-start_time)
        run_one_binary(mass)


run_grid()
