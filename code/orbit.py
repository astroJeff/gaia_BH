import numpy as np
from scipy.optimize import newton
from _skysub import parellipse


# Constants
Msun = 1.989e33  # in g
G = 6.674e-8  # in cgs
secday = 3600.0*24.0
secyer = 3600.0*24.0*365.25
AUincm = 1.496e13
pcincm = 3.086e18


def P_to_a(M1, M2, P):
    """ Orbital period (seconds) to separation (AU) """
    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    return np.power(mu/(n*n), 1.0/3.0) / AUincm

def a_to_P(M1, M2, a):
    """ Orbital separation (AU) to period (seconds) """
    mu = G * (M1 + M2)
    n = np.sqrt(mu/(a**3 * AUincm**3))
    return 2.0*np.pi / n

def get_f(p, t_set):
    """
    Calculate the true anomaly

    Arguments
    ---------
    p : tuple
        Set of model parameters
    t_set : float, numpy array
        Observations times for true anomaly calculation

    Returns
    -------
    f_set : float, numpy array
        Set of true anomalies
    """

    # Load model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate orbital separation
    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)

    # Kepler's equation
    def func_E(x,n,t,tau,e):
        return n*(t-tau) - x + e*np.sin(x)

    # Iterate through t_set, calculate eccentric anomaly - requires a root finder, so cannot be vectorized
    E_set = np.array([])
    if isinstance(t_set, np.ndarray):
        for t in t_set:
            try:
                # Use Newton root-finder to solve Kepler's equation
                E_new = newton(func_E, t, args=(n,t,tau,e))
            except:
                return None
            E_set = np.append(E_set, E_new)
    else:
        try:
            # Use Newton root-finder to solve Kepler's equation
            E_new = newton(func_E, t_set, args=(n,t_set,tau,e))
        except:
            return None
        E_set = np.append(E_set, E_new)

    # For each eccentric anomaly, calculate the true anomaly
    f_set = np.array([])
    for E in E_set:
        f = np.arccos((np.cos(E)-e)/(1.0-e*np.cos(E)))
        if np.sin(E) < 0:
            f = 2.0*np.pi - f

        f_set = np.append(f_set, f)

    return f_set



def get_r_dot(p, t):
    """
    Function to get the r component of the velocity

    Arguments
    ---------
    p : tuple
        Set of model parameters
    t : float, numpy array
        Observation times

    Returns
    -------
    r_dot : float, numpy array
        Radial velocities
    """

    # Load up model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate the orbital separation
    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)

    # Calculate true anomalies
    f = get_f(p, t)
    if f is None: return None

    # Calculate radial velocities
    r_dot = A * e * np.sin(f) / np.sqrt(1.0 - e*e) * n

    return r_dot


def get_r_f_dot(p, t):
    """
    Calculate the phi component of the velocity

    Arguments
    ---------
    p : tuple
        Set of model parameters
    t : float, numpy array
        Observation times

    Returns
    -------
    r_f_dot : float, numpy array
        Azimuthal velocities
    """

    # Load model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate the orbital separation
    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)

    # Calculate true anomalies
    f = get_f(p, t)
    if f is None: return None

    # Calculate the azimuthal velocities
    r_f_dot = n * A / np.sqrt(1.0 - e*e) * (1.0 + e*np.cos(f))

    return r_f_dot


def get_RV(p, t):
    """
    Calculate the radial velocities

    Arguments
    ---------
    p : tuple
        Set of model parameters
    t : float, numpy array
        Observation times

    Returns
    -------
    RV : float
        Radial Velocities
    """

    # Load model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate A
    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)

    # Get true anomalies
    f = get_f(p, t)
    if f is None: return None

    # Get two RV components (these are in cm/s)
    rv = n * A / np.sqrt(1.0 - e*e) * (np.cos(omega+f) * np.sin(I) + e * np.cos(omega) * np.sin(I))
    rv = rv * M2 / (M1+M2)

    # Convert from cm/s to km/s
    rv = rv / 1.0e5

    # Adjust for system's radial velocity
    rv = rv + gamma

    return rv




def get_ra_dec(p, t):
    """
    Function to get astrometry for a particular orbital solution and time

    Arguments
    ---------
    p : tuple
        Set of model parameters
    t : float, numpy array
        Observation times

    Returns
    -------
    ra, dec : output coordinates (deg)
    """

    # load model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    # Calculate orbital separation
    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)

    # Calculate true anomalies
    f = get_f(p, t)
    if f is None: return None, None

    # Get the distance from M2 to the barycenter
    sep = A * (1.0 - e**2) / (1.0 + e * np.cos(f)) * (M2 / (M1+M2))

    # Calculate the orbital motion
    d_alpha = sep * (np.cos(Omega)*np.cos(omega+f) - np.sin(Omega)*np.sin(omega+f)*np.cos(I)) / (distance*pcincm)
    d_delta = sep * (np.sin(Omega)*np.cos(omega+f) + np.cos(Omega)*np.sin(omega+f)*np.cos(I)) / (distance*pcincm)

    # Convert to degrees, add in center of motion position
    alpha_out = sys_ra + d_alpha * (180.0/np.pi)
    delta_out = sys_dec + d_delta * (180.0/np.pi)

    return alpha_out, delta_out





def get_ra_dec_all(p, t, include_orbit=True, include_parallax=True, include_proper_motion=True):
    """
    Obtain position coordinates of a star in a binary including orbit, parallax, and proper motion

    Arguments
    ---------
    p : tuple
        Set of orbital parameters
    t : float, numpy array
        Observation times

    Returns
    -------
    ra, dec : float, numpy array
        Output positions (deg)
    """

    # Load up model parameters
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance, pm_ra, pm_dec = p
    p_orb = sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance

    # Get parallax from distance
    parallax = 1.0e3/distance

    # Calculate the effect on position from parallax
    if include_parallax:
        sys_ra_hour = sys_ra / 360.0 * 24.0  # convert from decimal degrees to decimal hours

        # Calculate parallax factors using parellipse from  J. Thorstensen's code SKYCALC
        output = parellipse(t, sys_ra_hour, sys_dec, 2015.0, 0.0, 0.0)

        # Calculate actual effect on position from parallax by multiplying by parallax
        # Since parallax is in units of mas, positions are currently also in mas
        ra_out = output[0] * parallax
        dec_out = output[1] * parallax
    else:
        ra_out = np.zeros(len(t))
        dec_out = np.zeros(len(t))


    # Add the effect from proper motion (also in mas)
    if include_proper_motion:
        ra_out += pm_ra * (t*secday - tau)/secday/365.25
        dec_out += pm_dec * (t*secday - tau)/secday/365.25


    # Calculate the effects of orbital motion
    if include_orbit:
        ra_orb, dec_orb = get_ra_dec(p_orb, t*secday)
        if ra_orb is None or dec_orb is None: return None, None
    else:
        ra_orb, dec_orb = 0.0, 0.0

    # Combine the effects of parallax, proper motion, and orbital motion (deg)
    ra_out = ra_out/(3600.0*1.0e3) + ra_orb
    dec_out = dec_out/(3600.0*1.0e3) + dec_orb

    return ra_out, dec_out
