import numpy as np
import time
from scipy.interpolate import interp1d
from scipy.stats import multivariate_normal
import sys

import emcee
import pygaia.errors.photometric as gaia_photometric

import orbit
import prob
import gaia_tools as gaia
import const as c
import photometry





def load_t_obs(filename="../data/cadence/J1102.csv"):
    data = np.genfromtxt(filename, delimiter=',', names=True)

    return data




def start_position(M1, M2, distance, P_orb, t_obs, pm_ra, pm_dec,
                   sys_ra=165.5728333, sys_dec=41.2209444, Omega=np.pi/3.0,
                   omega=np.pi/3.0, I=np.pi/3.0, e=0.01, gamma=0.0):

    # Observation times
    tau = np.min(t_obs) * c.secday

    p = sys_ra, sys_dec, Omega, omega, I, tau, e, P_orb, gamma, M1, M2, distance, pm_ra, pm_dec

    return p



def run_one_M2(P_orb=100.0, M2=1.0e-3, nburn=1000, nrun=5000):

    dist_set = np.logspace(1, 4, 4)
    M1_set = [1.0, 10.0]
    # dist_set = [2.0]
    # M1_set = [1.0]

    for dist in dist_set:
        for M1 in M1_set:

            sampler = run_one_binary(dist, M1*c.Msun, M2*c.Msun, P_orb*c.secday, nburn=nburn, nrun=nrun)

            fileout = "../data/M1_" + str(int(M1)) + '_M2_%.3f'%M2 + '_dist_' + str(int(dist)) + '_Porb_%.3f'%P_orb
            np.save(fileout + "_chains.npy", sampler.chain[:,::10,:])
            np.save(fileout + "_acceptancefractions.npy", sampler.acceptance_fraction)



def run_one_binary(distance, M1, M2, P_orb, nburn=1000, nrun=5000, include_M2_photo=True):

    # Load Gaia observation times
    data = load_t_obs()
    t_obs = data['ObservationTimeAtBarycentreBarycentricJulianDateInTCB']
    obs_angle = data['scanAnglerad']

    # Set proper motion
    pm_ra = 10.0  # in mas/yr
    pm_dec = 10.0  # in mas/yr

    # Create tuple of initial values
    p = start_position(M1, M2, distance, P_orb, t_obs, pm_ra, pm_dec,
                       sys_ra=165.5728333, sys_dec=41.2209444, Omega=np.pi/3.0,
                       omega=np.pi/3.0, I=np.pi/3.0, e=0.01, gamma=0.0)

    # Get observation uncertainty - use Main Sequence star tables
    # get_V_IC, get_M_G = load_stellar_photometry()
    # V_IC = get_V_IC(M1/c.Msun)
    # G_mag = get_M_G(M1/c.Msun) + 5.0*np.log10(distance) - 5.0
    G_mag_obs = None
    G_mag_err = None
    V_IC = photometry.get_V_IC(M1/c.Msun)
    G_mag_obs = photometry.get_G_mag_apparent(M1/c.Msun, distance)
    G_mag_err = gaia_photometric.gMagnitudeErrorEoM(G_mag_obs, nobs=len(t_obs))

    AL_err = gaia.get_single_obs_pos_err(G=G_mag_obs, V_IC=V_IC, RA=165.5728333, Dec=41.2209444, DIST=distance)
    AL_err *= 1.0e3  # in mas, not asec

    if not include_M2_photo: G_mag_obs = None

    # Calculate observed positions
    ra_obs, dec_obs = get_pos_obs(p, t_obs, obs_angle, AL_err)

    # Initialize walkers using start position
    sampler, p0 = initialize_walkers(p, ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err, nwalkers=128)

    # Run sampler
    sampler = run_emcee(sampler, p0, nburn=nburn, nrun=nrun)

    return sampler



def initialize_walkers(p, ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err, nwalkers=128):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance, pm_ra, pm_dec = p

    # Define sampler
    sampler = emcee.EnsembleSampler(nwalkers=nwalkers,
                                    dim=13,
                                    lnpostfn=prob.get_ln_posterior,
                                    args=(ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err))

    # Initialize walkers
    sys_ra_set = sys_ra + 1.0e-10*np.random.normal(size=nwalkers)
    sys_dec_set = sys_dec + 1.0e-10*np.random.normal(size=nwalkers)
    Omega_set = Omega + 1.0e-5*np.random.normal(size=nwalkers)
    omega_set = omega + 1.0e-5*np.random.normal(size=nwalkers)
    I_set = I + 1.0e-5*np.random.normal(size=nwalkers)
    tau_set = tau + 1.0e3*np.random.normal(size=nwalkers)
    e_set = e + 1.0e-6*np.random.normal(size=nwalkers)
    P_set = P * (1. + 1e-3*np.random.normal(size=nwalkers))
    M1_set = M1 * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    M2_set = M2 * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    distance_set = distance * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    pm_ra_set = pm_ra * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    pm_dec_set = pm_dec * (1. + 1.0e-6*np.random.normal(size=nwalkers))

    p0 = np.array([sys_ra_set,
                   sys_dec_set,
                   Omega_set,
                   omega_set,
                   I_set,
                   tau_set,
                   e_set,
                   P_set,
                   M1_set,
                   M2_set,
                   distance_set,
                   pm_ra_set,
                   pm_dec_set]).T

    return sampler, p0



def get_pos_obs(p, t_obs, obs_angle, AL_err):

    ra_obs = np.zeros(len(t_obs))
    dec_obs = np.zeros(len(t_obs))

    for i, t in enumerate(t_obs):
        ra_obs_i, dec_obs_i = orbit.get_ra_dec_all(p, t)
        if ra_obs_i is None or dec_obs_i is None: return None, None

        ra_obs[i], dec_obs[i] = ra_obs_i, dec_obs_i

        # Covariance matrix
        obs_cov = np.array([[AL_err**2, 0.0],
                            [0.0, gaia.AC_err**2]])
        rotation_matrix = np.array([[np.cos(obs_angle[i]), -np.sin(obs_angle[i])],
                                    [np.sin(obs_angle[i]), np.cos(obs_angle[i])]])
        new_cov = rotation_matrix @ obs_cov @ rotation_matrix.T

        # Add multivariate errors
        obs_astrometry = multivariate_normal(mean=np.zeros(2), cov=new_cov)
        obs_sample = obs_astrometry.rvs(size=1)
        ra_obs[i] += obs_sample[0]/3600.0/1.0e3
        dec_obs[i] += obs_sample[1]/3600.0/1.0e3

    # Errors are in mas, while ra_tmp is in degrees
    # ra_obs = ra_tmp + (pos_err/3600.0/1.0e3)*np.random.normal(size=len(t_obs))
    # dec_obs = dec_tmp + (pos_err/3600.0/1.0e3)*np.random.normal(size=len(t_obs))

    return ra_obs, dec_obs


def run_emcee(sampler, p0, nburn=100, nrun=1000):

    pos,prob,state = sampler.run_mcmc(p0, N=nburn)

    pos,prob,state = sampler.run_mcmc(pos, N=nrun)

    return sampler





def run_only_one(dist, M1, M2, P_orb, nburn, nrun, include_M2_photo=False):

    sampler = run_one_binary(dist, M1*c.Msun, M2*c.Msun, P_orb*c.secday,
                             nburn=nburn, nrun=nrun, include_M2_photo=include_M2_photo)

    fileout = "../data/M1_" + str(int(M1)) + '_M2_%.3f'%M2 + '_dist_' + str(int(dist)) + '_Porb_%.3f'%P_orb
    if include_M2_photo: fileout = fileout + '_photo' 
    np.save(fileout + "_chains.npy", sampler.chain[:,::100,:])
    np.save(fileout + "_acceptancefractions.npy", sampler.acceptance_fraction)


# Binary parameters
dist = float(sys.argv[1])
M1 = float(sys.argv[2])
M2 = float(sys.argv[3])
P_orb = float(sys.argv[4])

include_M2_photo = False
if len(sys.argv) > 5: include_M2_photo=True

# Run the binary
# run_only_one(dist, M1, M2, P_orb, 1000, 10000, include_M2_photo=False)
run_only_one(dist, M1, M2, P_orb, 1000, 10000, include_M2_photo=include_M2_photo)








##############
#   OLD CODE #
##############
#
#
#
# def load_data(file_name):
#
#     data = np.genfromtxt(file_name)
#
#     names = ["time","ra","dec","ra_err","dec_err","rv","rv_err","plx","plx_err"]
#     obs_pos = np.recarray(len(data[:,0]), names=names,
#                           formats=['float64,float64,float64,float64,float64,float64,float64,float64,float64'])
#
#     obs_pos['time'] = data[:,0]  # in seconds
#     obs_pos['ra'] = data[:,1]/3600.0/1.0e3  # in degrees
#     obs_pos['dec'] = data[:,2]/3600.0/1.0e3  # in degrees
#     obs_pos['rv'] = (2./3.) * data[:,3]/1.0e3  # in km/s
#
#     obs_pos['ra_err'] = data[:,4]/1.0e6 # in asec
#     obs_pos['dec_err'] = data[:,4]/1.0e6 # in asec
#     obs_pos['rv_err'] = 1.0
#
#     return obs_pos
#
#
#
#
# def adjust_distance(obs_pos, distance_new, distance_current=2003.0, mag_current=12.5):
#     """Adjust object according to new distance
#
#     """
#
#     V_IC = -0.5417
#
#     mag_new = gaia.calc_new_magnitude(mag_current, distance_current, distance_new)
#
#     plx_err = gaia.get_plx_err(G=mag_new, V_IC=V_IC) * 1.0e3  # plx error in mas
#
#     obs_pos['plx'] = 1.0e3 / distance_new
#     obs_pos['plx_err'] = plx_err
#
#
#     # Adjust position coordinates
#     obs_pos['ra'] *= distance_current / distance_new
#     obs_pos['dec'] *= distance_current / distance_new
#
#     obs_pos['ra_err'] = gaia.get_single_obs_err(G=mag_new, V_IC=V_IC)
#     obs_pos['dec_err'] = gaia.get_single_obs_err(G=mag_new, V_IC=V_IC)
#
#     return obs_pos
#
# def add_uncertainties(obs_pos):
#     obs_pos['ra'] += np.random.normal(scale=obs_pos['ra_err']/3600.0, size=len(obs_pos))
#     obs_pos['dec'] += np.random.normal(scale=obs_pos['dec_err']/3600.0, size=len(obs_pos))
#     obs_pos['rv'] += np.random.normal(scale=obs_pos['rv_err'], size=len(obs_pos))
#
#     return obs_pos
#
#
# def get_truths(distance_new):
#
#     sys_ra = 1.0e-10
#     sys_dec = 1.0e-10
#     Omega = 2.1
#     omega = 2.5
#     I = 0.3
#     # tau = 0.33 + P/2.1
#     e = 0.7
#     P = 2.0*3600.0*24.0*365.25
#     tau = P/500.0
#     tau = 100.0
#     gamma = 1.0e-4
#     M1 = 10.0*c.Msun
#     M2 = 20.0*c.Msun
#
#     truths = sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance_new
#
#     return truths
#
#
# def initialize_sampler(truths, obs_pos, nwalkers=100):
#
#     ndim = 12
#
#     # Starting array
#     p0 = np.zeros((nwalkers, ndim))
#
#     # Random deviations from best fit
#     delta = 1.0 + np.random.normal(0.0, 0.001, size=(nwalkers, ndim))
#
#     # Adjust starting array accordingly
#     for i in range(ndim): p0[:,i] = truths[i] * delta[:,i]
#
#     return p0
#
#
# def run_emcee(p0, obs_pos):
#
#     nwalkers, ndim = p0.shape
#
#     # Define sampler
#     sampler = emcee.EnsembleSampler(nwalkers=nwalkers, dim=ndim, lnpostfn=prob.ln_posterior, args=(obs_pos,))
#
#     # Run burn-in
#     pos,ln_prob,state = sampler.run_mcmc(p0, N=1000)
#
#     # Run sampler
#     sampler.reset()
#     pos,ln_prob,state = sampler.run_mcmc(pos, N=5000)
#
#     # Downsample by factor of 10
#     chains = sampler.chain[:,::10,:]
#
#     return chains
#
#
# def run_one_binary(distance, rv_err):
#
#     # Load binary data
#     obs_pos = load_data("../data/Eccentric_True_RA_Dec_RV.txt")
#
#     # Adjust binary to new distance
#     obs_pos = adjust_distance(obs_pos, distance)
#
#     # Adjust binary's RV error
#     obs_pos['rv_err'] = rv_err
#
#     # Add uncertainties to position and RV measurements
#     obs_pos = add_uncertainties(obs_pos)
#
#     # Calculate truths for the binary
#     truths = get_truths(distance)
#
#     # Initialize the emcee sampler
#     p0 = initialize_sampler(truths, obs_pos, nwalkers=100)
#
#     # Run the sampler
#     chains = run_emcee(p0, obs_pos)
#
#     # Save the chains
#     file_out = "../data/binary_chains/distance" + "%.2f" % distance + "_rverr" + "%.2f" % rv_err + ".npy"
#     np.save(file_out, chains)
#
#
#
# def run_grid(n_dist=10, n_rv_err=10, dist_range=(50, 10000), rv_err_range=(0.3, 20.0)):
#
#     start_time = time.time()
#
#     dist_set = 10**np.linspace(np.log10(dist_range[0]), np.log10(dist_range[1]), n_dist)
#     rv_err_set = 10**np.linspace(np.log10(rv_err_range[0]), np.log10(rv_err_range[1]), n_rv_err)
#
#     for dist in dist_set:
#         for rv_err in rv_err_set:
#
#             print("Running:", dist, rv_err, time.time()-start_time)
#             run_one_binary(dist, rv_err)
#
#
# run_grid()
