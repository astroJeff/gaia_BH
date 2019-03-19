import numpy as np
from scipy.stats import multivariate_normal

import orbit
import const as c
from gaia_tools import AC_err

import photometry



def get_ln_likelihood(p, ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err):
    """ Calculate the natural log of the likelihood for a set of model parameters

    Arguments
    ---------
    p : tuple
        Set of model parameters
    ra_obs : nd_array
        Set of positions in right ascension (deg)
    dec_obs : nd_array
        Set of positions in declination (deg)
    t_obs : nd_array
        Set of observation times (Julian Days)
    obs_angle : nd_array
        Set of position angles that Gaia rotates on the sky (rad)
    AL_err : float
        Along-scan error (mas)
    G_mag_obs : float
        Gaia G band magnitude
    G_mag_err : float
        Gaia G band magnitude error

    Returns
    -------
    ln_likelihood : float
        Natural log of the likelihood

    """

    # Load the set of model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance, pm_ra, pm_dec = p

    ln_likelihood = 0.0

    # Create zero arrays
    ra_tmp = np.zeros(len(t_obs))
    dec_tmp = np.zeros(len(t_obs))

    # Iterate through individual observations
    for i, t in enumerate(t_obs):

        # Calculate the orbital position as a function of model parameters at a time t
        ra_tmp_i, dec_tmp_i = orbit.get_ra_dec_all(p, t)
        if ra_tmp_i is None or dec_tmp_i is None: return -np.inf
        ra_tmp[i], dec_tmp[i] = ra_tmp_i, dec_tmp_i

        # Create the array of observed positions
        obs_pos_arr = np.array([ra_obs[i], dec_obs[i]])

        # Create the array of calculated positions
        pos_arr = np.array([ra_tmp[i], dec_tmp[i]])

        # Create the along-scan vs. across-scan covariance matrix (in deg)
        obs_cov = np.array([[(AL_err/1.0e3/3600.0)**2, 0.0],
                            [0.0, (AC_err/1.0e3/3600.0)**2]])

        # Calculate the rotation matrix corresponding to the angle observed by Gaia
        rotation_matrix = np.array([[np.cos(obs_angle[i]), -np.sin(obs_angle[i])],
                                    [np.sin(obs_angle[i]), np.cos(obs_angle[i])]])

        # Obtain the covariance matrix in ra vs. dec space
        new_cov = rotation_matrix @ obs_cov @ rotation_matrix.T

        # Create the position object for the observation
        obs_astrometry = multivariate_normal(mean=pos_arr, cov=new_cov)

        # Calculate the log-likelihood for uncertainties in ra-dec space
        ln_likelihood += np.log(obs_astrometry.pdf(obs_pos_arr))

    # Include limits on the photometry
    if G_mag_obs is not None:

        # Obtain the apparent magnitude of the modeled star given its mass and distance
        G_mag_model = photometry.get_G_mag_apparent(M1/c.Msun, distance)

        # Adjust the log of the likelihood based on photometry
        ln_likelihood += -(G_mag_obs - G_mag_model)**2/(2.0*G_mag_err**2)

    return ln_likelihood


def get_ln_prior(p):
    """ Calculate the prior probability of a set of model parameters

    Arguments
    ---------
    p : tuple
        Model parameters

    Returns
    -------
    ln_prior : float
        Natural log of the prior probabilty for the set of model parameters
    """

    # Load the set of model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, M1, M2, distance, pm_ra, pm_dec = p

    ln_prior = 0.0

    # Limit the range of model parameters
    if sys_ra < 0.0 or sys_ra > 360.0: return -np.inf
    if sys_dec < -90.0 or sys_dec > 90.0: return -np.inf
    if Omega < 0.0 or Omega > 2.0*np.pi: return -np.inf
    if omega < 0.0 or omega > np.pi: return -np.inf
    if I < 0.0 or I > np.pi: return -np.inf
    if tau < 2451545.*c.secday or tau > 2471545.*c.secday: return -np.inf
    if e < 0. or e >= 1.0: return -np.inf
    if P < 0.0 or P > 2.0e4*c.secday: return -np.inf
    if M1 < 0.0 or M1 > 20.0*c.Msun: return -np.inf
    if M2 < 0.0 or M2 > 20.0*c.Msun: return -np.inf
    if distance < 0.0 or distance > 20000.0: return -np.inf
    if pm_ra < -1.0e4 or pm_ra > 1.0e4: return -np.inf
    if pm_dec < -1.0e4 or pm_dec > 1.0e4: return -np.inf

    # Log flat priors on mass
    ln_prior -= np.log(M1/c.Msun)
    ln_prior -= np.log(M2/c.Msun)

    # Prior on inclination angle
    ln_prior += np.log(np.sin(I)/2.0)

    # Distance prior Lutz-Kelker bias
    distance_max = 20.0e3  # Maximum distance of 20 kpc
    ln_prior += np.log(3.0 * distance**2 / distance_max**3)

    return ln_prior


def get_ln_posterior(p, ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err):
    """ Calculate the posterior probability of a set of model parameters

    Arguments
    ---------
    p : tuple
        Set of model parameters
    ra_obs : nd_array
        Set of positions in right ascension (deg)
    dec_obs : nd_array
        Set of positions in declination (deg)
    t_obs : nd_array
        Set of observation times (Julian Days)
    obs_angle : nd_array
        Set of position angles that Gaia rotates on the sky (rad)
    AL_err : float
        Along-scan error (mas)
    G_mag_obs : float
        Gaia G band magnitude
    G_mag_err : float
        Gaia G band magnitude error

    Returns
    -------
    ln_posterior : float
        Natural log of the posterior probability
    """

    # Load set of model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, M1, M2, distance, pm_ra, pm_dec = p

    # Natural log of the prior
    ln_prior = get_ln_prior(p)
    if np.isinf(lp): return -np.inf

    # We don't model the radial velocity zero point, but its required by the functions
    gamma = 0.0 # Holder variable - not a useful parameter
    p = sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance, pm_ra, pm_dec

    # Natural log of the likelihood
    ln_likelihood = get_ln_likelihood(p, ra_obs, dec_obs, t_obs, obs_angle, AL_err, G_mag_obs, G_mag_err)
    if np.isnan(ll): return -np.inf

    # Return the natural log of the posterior function
    return ln_prior + ln_likelihood
