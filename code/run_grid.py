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
    """
    Load up Gaia cadence data

    Arguments
    ---------
    filename : string
        File path and name for gaia cadence data

    Returns
    -------
    data : numpy array
        Gaia observation cadence data
    """

    data = np.genfromtxt(filename, delimiter=',', names=True)

    return data




def start_position(M1, M2, distance, P_orb, t_obs, pm_ra, pm_dec,
                   sys_ra=165.5728333, sys_dec=41.2209444, Omega=np.pi/3.0,
                   omega=np.pi/3.0, I=np.pi/3.0, e=0.01, gamma=0.0):
    """
    Set the starting position for model parameters

    Arguments
    ---------
    M1, M2 : float
        Component masses (Msun)
    distance : float
        Distance to the binary (pc)
    P_orb : float
        Orbital period (days)
    t_obs : float
        Set of observation times
    pm_ra, pm_dec : float
        Proper motions (mas/yr)
    sys_ra, sys_dec : float
        Center of mass position coordinates
    Omega : float
        Longitude of the ascending node
    omega : float
        Argument of periapse
    I : float
        Inclination angle
    e : float
        Eccentricity
    gamma : float
        radial velocity zero point

    Returns
    -------
    p : tuple
        Set of model parameters
    """

    # Observation times
    tau = np.min(t_obs) * c.secday

    p = sys_ra, sys_dec, Omega, omega, I, tau, e, P_orb, gamma, M1, M2, distance, pm_ra, pm_dec

    return p



def run_one_M2(P_orb=100.0, M2=1.0e-3, nburn=1000, nrun=5000):
    """
    Run grid of models in distance and M1. Save model results.

    Parameters
    ----------
    P_orb : float
        Orbital period (days)
    M2 : float
        Companion mass (Msun)
    nburn : int
        Number of burn-in steps
    nrun : int
        Number of steps to run
    """

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



def run_one_binary(distance, M1, M2, P_orb, nburn=1000, nrun=5000, include_M2_photo=True, include_RV=False,
                   t_obs_filename="../data/cadence/J1102.csv", fileout=None):
    """
    Run sampler on one binary

    Arguments
    ---------
    distance : float
        Distance to the binary
    M1, M2 : float
        Component masses (Msun)
    P_orb : float
        Orbital period (days)
    nburn, nrun : int
        Number of steps in burn-in and run
    include_M2_photo : bool
        Add photometric constraints to likelihood
    include_RV : bool
        Add radial velocity constraints to likelihood

    Returns
    -------
    sampler : emcee sampler
        MCMC emcee sampler
    """

    # Load Gaia observation times and angles
    data = load_t_obs(filename=t_obs_filename)
    t_obs = data['ObservationTimeAtBarycentreBarycentricJulianDateInTCB']
    obs_angle = data['scanAnglerad']

    # Set proper motion
    pm_ra = 10.0  # in mas/yr
    pm_dec = 10.0  # in mas/yr

    # Create tuple of initial values
    p = start_position(M1, M2, distance, P_orb, t_obs, pm_ra, pm_dec,
                       sys_ra=165.5728333, sys_dec=41.2209444, Omega=np.pi/3.0,
                       omega=np.pi/3.0, I=np.pi/3.0, e=0.01, gamma=0.0)

    # Gaia G mmagnitude and error
    G_mag_obs = None
    G_mag_err = None

    # Obtain colors, Gaia G magnitude and error from tables
    # V_IC = photometry.get_V_IC(M1/c.Msun)
    G_mag_obs = photometry.get_G_mag_apparent(M1/c.Msun, distance)
    G_mag_err = gaia_photometric.gMagnitudeErrorEoM(G_mag_obs, nobs=len(t_obs))

    # Calculate Gaia along-scan pointing precision
    AL_err = gaia.get_single_obs_pos_err_limit(g_mag=G_mag_obs)
    # AL_err = gaia.get_single_obs_pos_err(G=G_mag_obs, V_IC=V_IC, RA=165.5728333, Dec=41.2209444, DIST=distance)
    AL_err *= 1.0e3  # in mas, not asec

    # Include photometric constraints?
    if not include_M2_photo: G_mag_obs = None

    # Include RV constraints
    if not include_RV:
        t_obs_RV = None
        RV_obs = None
        RV_err = None
    else:
        t_obs_RV = data['ObservationTimeAtBarycentreBarycentricJulianDateInTCB'][data['CcdRow17'] <= 4] # RV obs only on CCDs 1-4
        RV_obs, RV_err = get_RV_obs(p, t_obs_RV, RV_err=1.0)

    # Calculate observed positions
    ra_obs, dec_obs = get_pos_obs(p, t_obs, obs_angle, AL_err)

    # Initialize walkers using start position
    sampler, p0 = initialize_walkers(p, ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err,
                                     nwalkers=128, include_RV=include_RV, t_obs_RV=t_obs_RV, RV_obs=RV_obs, RV_err=RV_err)

    # Run sampler in batches of 1000
    batch = 1000
    for i in range(int(nrun/batch)+1):

        print(i, fileout)

        sampler, p0 = run_emcee(sampler, p0, nrun=np.min([batch, nrun - i*batch]))

        # Save sampler after every batch
        if fileout is not None:
            np.save(fileout + "_chains.npy", sampler.chain[:,::100,:])
            np.save(fileout + "_acceptancefractions.npy", sampler.acceptance_fraction)

    return sampler



def initialize_walkers(p, ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err,
                       nwalkers=128, include_RV=False, t_obs_RV=None, RV_obs=None, RV_err=None):
    """
    Initialize walkers in an n-ball around set of model parameters and emcee sampler

    Arguments
    ---------
    p : tuple
        Set of model parameters
    ra_obs, dec_obs : numpy array
        numpy array of positions
    t_obs : numpy array
        numpy array of observation times
    obs_angle : numpy array
        position angles observed by Gaia
    AL_err : float
        Along-scan single-pointing precision
    G_mag_obs, G_mag_err : float
        Gaia G magnitude and error
    nwalkers : int
        Number of walkers for emcee

    Parameters
    ----------
    sampler : emcee sampler
        Initialized emcee sampler
    p0 : numpy array
        Set of numpy arrays
    """

    # Load model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance, pm_ra, pm_dec = p

    # Vary probabilities depending radial velocities
    if include_RV:
        dim = 14
        func_posterior = prob.get_ln_posterior_RV
        args = (ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err, t_obs_RV, RV_obs, RV_err)
    else:
        dim = 13
        func_posterior = prob.get_ln_posterior
        args = (ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err)

    # Define sampler
    sampler = emcee.EnsembleSampler(nwalkers=nwalkers,
                                    dim=dim,
                                    lnpostfn=func_posterior,
                                    args=args)

    # Initialize walkers
    sys_ra_set = sys_ra + 1.0e-10*np.random.normal(size=nwalkers)
    sys_dec_set = sys_dec + 1.0e-10*np.random.normal(size=nwalkers)
    Omega_set = Omega + 1.0e-5*np.random.normal(size=nwalkers)
    omega_set = omega + 1.0e-5*np.random.normal(size=nwalkers)
    I_set = I + 1.0e-5*np.random.normal(size=nwalkers)
    tau_set = tau + 1.0e3*np.random.normal(size=nwalkers)
    e_set = e + 1.0e-6*np.random.normal(size=nwalkers)
    P_set = P * (1. + 1e-3*np.random.normal(size=nwalkers))
    if include_RV: gamma_set = gamma + 1.0e-2*np.random.normal(size=nwalkers)
    M1_set = M1 * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    M2_set = M2 * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    distance_set = distance * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    pm_ra_set = pm_ra * (1. + 1.0e-6*np.random.normal(size=nwalkers))
    pm_dec_set = pm_dec * (1. + 1.0e-6*np.random.normal(size=nwalkers))

    # Numpy array of walkers
    if include_RV:
        p0 = np.array([sys_ra_set,
                       sys_dec_set,
                       Omega_set,
                       omega_set,
                       I_set,
                       tau_set,
                       e_set,
                       P_set,
                       gamma_set,
                       M1_set,
                       M2_set,
                       distance_set,
                       pm_ra_set,
                       pm_dec_set]).T
    else:
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
    """
    Generate mock observations

    Arguments
    ---------
    p : tuple
        Set of model parameters
    t_obs : numpy array
        Observation times
    obs_angles : numpy array
        Position angles of observations (rad)
    AL_err : float
        Along-scan single-pointing precision

    Returns
    -------
    ra_obs, dec_obs : float
        Observed positions (deg)
    """

    # RA and Dec positions
    ra_obs = np.zeros(len(t_obs))
    dec_obs = np.zeros(len(t_obs))

    # Iterate through observation times
    for i, t in enumerate(t_obs):

        # Calculate the orbital position as a function of model parameters at a time t
        ra_obs_i, dec_obs_i = orbit.get_ra_dec_all(p, t)
        if ra_obs_i is None or dec_obs_i is None: return None, None

        # Create the array of observed positions
        ra_obs[i], dec_obs[i] = ra_obs_i, dec_obs_i

        # Create the along-scan vs. across-scan covariance matrix (in deg)
        obs_cov = np.array([[AL_err**2, 0.0],
                            [0.0, gaia.AC_err**2]])

        # Calculate the rotation matrix corresponding to the angle observed by Gaia
        rotation_matrix = np.array([[np.cos(obs_angle[i]), -np.sin(obs_angle[i])],
                                    [np.sin(obs_angle[i]), np.cos(obs_angle[i])]])

        # Obtain the covariance matrix in ra vs. dec space
        new_cov = rotation_matrix @ obs_cov @ rotation_matrix.T

        # Create the position object for the observation
        obs_astrometry = multivariate_normal(mean=np.zeros(2), cov=new_cov)

        # Add multivariate errors
        obs_sample = obs_astrometry.rvs(size=1)
        ra_obs[i] += obs_sample[0]/3600.0/1.0e3
        dec_obs[i] += obs_sample[1]/3600.0/1.0e3

    return ra_obs, dec_obs


def get_RV_obs(p, t_obs_RV, RV_err=1.0):
    """
    Get RV samples

    Arguments
    ---------
    p : tuple
        Set of model parameters
    t_obs_RV : nd_array
        Set of radial velocity observation times
    RV_err : float
        Radial velocity error (km/s)

    Returns
    -------
    RV_obs : nd_array
        Array of radial velocities
    RV_err : float
        Radial velocity error
    """

    sys_ra, sys_dec, Omega, omega, I, tau, e, P_orb, gamma, M1, M2, distance, pm_ra, pm_dec = p
    p_RV = sys_ra, sys_dec, Omega, omega, I, tau, e, P_orb, gamma, M1, M2, distance

    # Calculate RVs
    RV_obs = orbit.get_RV(p_RV, t_obs_RV*c.secday)

    # Add measurement errors
    RV_obs = RV_obs + np.random.normal(0.0, RV_err, len(t_obs_RV))

    return RV_obs, RV_err


def run_emcee(sampler, p0, nrun=1000):
    """
    Run the emcee sampler

    Arguments
    ---------
    sampler : emcee sampler
        Initialized emcee sampler object
    p0 : numpy array
        Set of initialized walkers
    nrun : int
        Number of run steps

    Returns
    -------
    sampler : emcee sampler
        Run emcee sampler object
    """

    pos,prob,state = sampler.run_mcmc(p0, N=nrun)

    return sampler, prob





def run_only_one(dist, M1, M2, P_orb, nrun, include_M2_photo=False, include_RV=False, t_obs_filename="../data/cadence/J1102.csv"):
    """
    Run a single binary

    Arguments
    ---------
    dist : float
        Distance to the binary (pc)
    M1, M2 : float
        Component masses (Msun)
    P_porb : float
        Orbital period (days)
    nburn, nrun : int
        Number of burn-in and run steps
    include_M2_photo : bool
        Include photometric constraints
    include_RV : bool
        Include radial velocity constraints
    t_obs_filename : string
        Filename and path to observation time data
    """

    # Create the fileout string
    fileout = "../data/M1_" + str(int(M1)) + '_M2_%.3f'%M2 + '_dist_' + str(int(dist)) + '_Porb_%.3f'%P_orb
    if include_M2_photo: fileout = fileout + '_photo'
    if include_RV: fileout = fileout + '_RV'

    # Run a single binary
    sampler = run_one_binary(dist, M1*c.Msun, M2*c.Msun, P_orb*c.secday,
                             nburn=nburn, nrun=nrun, include_M2_photo=include_M2_photo,
                             include_RV=include_RV, t_obs_filename=t_obs_filename, fileout=fileout)

    # # Save sampler
    # np.save(fileout + "_chains.npy", sampler.chain[:,::100,:])
    # np.save(fileout + "_acceptancefractions.npy", sampler.acceptance_fraction)


# Binary parameters
dist = float(sys.argv[1])
M1 = float(sys.argv[2])
M2 = float(sys.argv[3])
P_orb = float(sys.argv[4])

include_M2_photo = False
if len(sys.argv) > 5:
    if float(sys.argv[5]) == 1.0: include_M2_photo = True

include_RV = False
if len(sys.argv) > 6:
    if float(sys.argv[6]) == 1.0: include_RV = True

filename = "../data/cadence/J1102.csv"
if len(sys.argv) > 7:
    filename = sys.argv[7]


# Run the binary
run_only_one(dist, M1, M2, P_orb, 40000, include_M2_photo=include_M2_photo, include_RV=include_RV, t_obs_filename=filename)
