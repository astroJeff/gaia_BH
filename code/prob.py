import numpy as np
import orbit


def ln_likelihood_rv(p, obs_rv):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    rv_model = orbit.get_RV(p, obs_rv["time"])

    # P( delta t_cool | M1_zams_t, M2_zams_t)
    rv_diff = rv_model - obs_rv["rv"]

    ln_P_temp = - rv_diff**2/(2.0*obs_rv["rv_err"]**2)

    return np.sum(ln_P_temp)

def ln_likelihood_ra(p, obs_pos):
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate position
    obs_time = obs_pos["time"]
    ra_out, dec_out = orbit.get_ra_dec(p, obs_time)

    # Work in degrees
    ra_diff = ra_out - obs_pos["ra"]

    ln_prob_ra = -np.log(obs_pos["ra_err"] * np.sqrt(2.0*np.pi)) - ra_diff**2/(2.0*obs_pos["ra_err"]**2)

    return np.sum(ln_prob_ra)

def ln_likelihood_dec(p, obs_pos):
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate position
    obs_time = obs_pos["time"]
    ra_out, dec_out = orbit.get_ra_dec(p, obs_time)

    # Work in asec-space
    dec_diff = dec_out - obs_pos["dec"]

    ln_prob_dec = -np.log(obs_pos["dec_err"] * np.sqrt(2.0*np.pi)) - dec_diff**2/(2.0*obs_pos["dec_err"]**2)

    return np.sum(ln_prob_dec)

def ln_likelihood_pos(p, obs_pos):
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate position
    obs_time = obs_pos["time"]
    ra_out, dec_out = orbit.get_ra_dec(p, obs_time)

    # Work in asec-space
    ra_diff = ra_out - obs_pos["ra"]
    dec_diff = dec_out - obs_pos["dec"]

    ln_prob_ra = - ra_diff**2/(2.0*obs_pos["ra_err"]**2)
    ln_prob_dec = - dec_diff**2/(2.0*obs_pos["dec_err"]**2)


    return np.sum(ln_prob_ra) + np.sum(ln_prob_dec)


def ln_likelihood_plx(p, obs_pos):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Only include first value here
    # Parallax in as, so distance in pc
    return - (obs_pos['plx'][0]-1.0/distance)**2/(2.0*obs_pos['plx_err'][0]**2)


def ln_likelihood_M2_photo(p, obs_pos, M2_obs, M2_obs_err):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    return - (M2_obs-M2)**2/(2.0*M2_obs_err**2)


# The posterior function calls the likelihood and prior functions
def ln_posterior(p, obs_pos, M2_obs, M2_obs_err):


    # Calculate prior on model first
    lp = ln_prior(p)

    # Calculate likelihood
    ll = ln_likelihood(p, obs_pos, M2_obs, M2_obs_err)

    return lp + ll



def ln_likelihood(p, obs_pos, M2_obs, M2_obs_err):

    ll = 0.0
    ll = ll + ln_likelihood_rv(p, obs_pos)
    ll = ll + ln_likelihood_pos(p, obs_pos)
    ll = ll + ln_likelihood_plx(p, obs_pos)
    ll = ll + ln_likelihood_M2_photo(p, obs_pos, M2_obs, M2_obs_err)

    if np.any(np.isnan(ll) | np.isinf(ll)):
        return -np.inf

    return ll




# Prior on parameters
def ln_prior(p):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    lp = 0.

    # Angles bound between 0 and 2pi
    if Omega < 0. or Omega > 2.0*np.pi:
        return -np.inf
    if omega < 0. or omega > 2.0*np.pi:
        return -np.inf

    # Needs to discern between left-handed and right-handed
    if I < 0. or I > np.pi:
        return -np.inf

    # Zero point needs to be constrained to first orbital period
    if tau < 0. or tau > P:
        return -np.inf

    # Eccentricity keeps the system bound
    if e < 0. or e > 1.0:
        return -np.inf

    # Period
    if P < 1.0e4 or P > 1.0e10:
        return -np.inf

    # Gamma
    if gamma < -1000.0 or gamma > 1000.0:
        return -np.inf

    # Constrain mass around 1 Msun
    if M1 < 0.1*orbit.Msun or M1 > 100.0*orbit.Msun:
        return -np.inf
    if M2 < 0.1*orbit.Msun or M2 > 100.0*orbit.Msun:
        return -np.inf

    # Prior on inclination angle
    lp = lp + np.log(np.sin(I)/2.0)

    # Prior on distance
    if distance < 1.0 or distance > 20000.0:
        return -np.inf

    if np.isnan(lp):
        return -np.inf

    return lp
