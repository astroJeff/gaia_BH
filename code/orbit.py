import numpy as np
from scipy.optimize import newton


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

# Function to get the true anomaly
def get_f(p, t_set):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)


    def func_E(x,n,t,tau,e):
        return n*(t-tau) - x + e*np.sin(x)

    E_set = np.array([])
    if isinstance(t_set, np.ndarray):
        for t in t_set:
            E_set = np.append(E_set, newton(func_E, t, args=(n,t,tau,e)))
    else:
        E_set = np.append(E_set, newton(func_E, t_set, args=(n,t_set,tau,e)))

    f_set = np.array([])
    for E in E_set:
        f = np.arccos((np.cos(E)-e)/(1.0-e*np.cos(E)))
        if np.sin(E) < 0:
            f = 2.0*np.pi - f

        f_set = np.append(f_set, f)


    return f_set


# Function to get the r component of the velocity
def get_r_dot(p, t):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)
    f = get_f(p, t)

    r_dot = A * e * np.sin(f) / np.sqrt(1.0 - e*e) * n

    return r_dot


# Function to get the phi component of the velocity (?)
def get_r_f_dot(p, t):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)
    f = get_f(p, t)

    r_f_dot = n * A / np.sqrt(1.0 - e*e) * (1.0 + e*np.cos(f))

    return r_f_dot

# Radial velocity
def get_RV(p, t):

    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)
    f = get_f(p, t)

    # r_dot = get_r_dot(p, t)
    # r_f_dot = get_r_f_dot(p, t)

    # Get two RV components (these are in cm/s)
    rv = n * A / np.sqrt(1.0 - e*e) * (np.cos(omega+f) * np.sin(I) + e * np.cos(omega) * np.sin(I))
    rv = rv * M2 / (M1+M2)

    # rv_1 = r_dot * M2 / (M1+M2) * (np.cos(Omega)*np.cos(omega+f) - np.sin(Omega)*np.sin(omega+f)*np.cos(I))
    # rv_2 = r_f_dot * M2 / (M1+M2) * (np.cos(Omega)*np.sin(omega+f) + np.sin(Omega)*np.cos(omega+f)*np.cos(I))

    # RV is in km/s
    # rv = -(rv_1 - rv_2) / 1.0e5
    rv = rv / 1.0e5

    # Adjust for system's radial velocity
    rv = rv + gamma

    return rv

# Function to get astrometry for a particular orbital solution and time

def get_ra_dec(p, t):
    sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance = p

    mu = G * (M1 + M2)
    n = 2.0*np.pi / P
    A = np.power(mu/(n*n), 1.0/3.0)
    f = get_f(p, t)

    # Get the distance from M2 to the barycenter
    sep = A * (1.0 - e**2) / (1.0 + e * np.cos(f)) * (M2 / (M1+M2))

    d_alpha = sep * (np.cos(Omega)*np.cos(omega+f) - np.sin(Omega)*np.sin(omega+f)*np.cos(I)) / (distance*pcincm)
    d_delta = sep * (np.sin(Omega)*np.cos(omega+f) + np.cos(Omega)*np.sin(omega+f)*np.cos(I)) / (distance*pcincm)

    # Work in arcseconds
    alpha_out = sys_ra + d_alpha * (180.0/np.pi)
    delta_out = sys_dec + d_delta * (180.0/np.pi)

    return alpha_out, delta_out
