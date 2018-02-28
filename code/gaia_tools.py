import numpy as np
from astropy.coordinates import SkyCoord
from astropy import units as u
from scipy.interpolate import interp1d
from pygaia.errors.astrometric import positionError as posErr
from pygaia.errors.astrometric import parallaxError as plxErr
from pygaia.errors.utils import averageNumberOfTransits as N_transit_ave


def get_single_obs_plx_err(G=None, V=None, V_IC=None, RA=None, Dec=None):
    """ A tool to determine the single observation astrometric precision. This is
    calculated in a bit of an ad hoc way. We call the overall parallax precision, then
    multiply by sqrt(N_obs), as this is the number that should have been used to determine
    the parallax precision. NOTE: this uses the Gaia DPAC package: pyGaia

    Inputs
    ------
    G or V : float (provide one of these)
    Gaia G-band or Visual V-band magnitude of the star

    V_IC : float
    Color index of the stellar type

    RA : float
    Right ascension in degrees

    Dec: float
    Declination in degrees

    Returns
    -------
    plx_err_single : float
    parallax uncertainty in arcseconds for a single observation

    """


    if V_IC is None:
        print("You must provide a color index in the form of V_IC.")
        return

    if G is None and V is None:
        print("You must provide at least one color: G or V")
        return
    
    if RA is None and Dec in None:
        print("You must provide a RA and Dec")
        return

    # KMB: get ecliptic latitude from ra/dec using astropy
    c = SkyCoord(ra=sys_ra*u.degree, dec=sys_dec*u.degree, frame='icrs')
    ecl_lat = c.GeocentricTrueEcliptic.lat

    # KMB: compute average N_obs for given ecliptic latitude
    N_obs_ave = N_transit_ave(ecl_lat)

    plx_err = plxErr(G=G, vmini=V_IC, beta=ecl_lat)*1.0e-6

    m = 1.2
    g_pi = 1.47/np.sin(ecl_lat)

    plx_err_single = np.sqrt(N_obs_ave) / (m * g_pi) * plx_err

    return plx_err_single



def get_single_obs_pos_err(G=None, V=None, V_IC=None, RA=None, Dec=None):
    """ A tool to determine the single observation astrometric precision. This is
    calculated in a bit of an ad hoc way. We call the overall parallax precision, then
    multiply by sqrt(N_obs), as this is the number that should have been used to determine
    the parallax precision. NOTE: this uses the Gaia DPAC package: pyGaia

    Inputs
    ------
    G or V : float (provide one of these)
    Gaia G-band or Visual V-band magnitude of the star

    V_IC : float
    Color index of the stellar type

    RA : float
    Right ascension in degrees

    Dec: float
    Declination in degrees    
     
    Returns
    -------
    pos_err_single : float
    position uncertainty in arcseconds for a single observation

    """

    if V_IC is None:
        print("You must provide a color index in the form of V_IC.")
        return

    if G is None and V is None:
        print("You must provide at least one color: G or V")
        return

    if RA is None and Dec in None:
        print("You must provide a RA and Dec")
        return

    # KMB: get ecliptic latitude from ra/dec using astropy
    c = SkyCoord(ra=sys_ra*u.degree, dec=sys_dec*u.degree, frame='icrs')
    ecl_lat = c.GeocentricTrueEcliptic.lat

    # KMB: compute average N_obs for given ecliptic latitude
    N_obs_ave = N_transit_ave(ecl_lat)

    pos_err = posErr(G=G, vmini=V_IC, beta=ecl_lat)*1.0e-6

    m = 1.2
    g_pi = 1.47/np.sin(ecl_lat)

    pos_err_single = np.sqrt(N_obs_ave) / (m * g_pi) * pos_err

    return pos_err_single


def get_plx_err(G=None, V=None, V_IC=None, RA=None, Dec=None):
    """ A tool to determine the parallax uncertainty for a star

    Inputs
    ------
    G or V : float (provide one of these)
    Gaia G-band or Visual V-band magnitude of the star

    V_IC : float
    Color index of the stellar type

    RA : float
    Right ascension in degrees

    Dec: float
    Declination in degrees

    Returns
    -------
    plx_err : float
    parallax uncertainty in arcseconds

    """


    if V_IC is None:
        print("You must provide a color index in the form of V_IC.")
        return

    if G is None and V is None:
        print("You must provide at least one color: G or V")
        return

    if G is None:  G = V - 0.0257 - 0.0924 * (V_IC) - 0.1623 * (V_IC)**2 + 0.0090 * (V_IC)**3

    # KMB: get ecliptic latitude from ra/dec using astropy
    c = SkyCoord(ra=sys_ra*u.degree, dec=sys_dec*u.degree, frame='icrs')
    ecl_lat = c.GeocentricTrueEcliptic.lat

    # KMB: compute average N_obs for given ecliptic latitude
    N_obs_ave = N_transit_ave(ecl_lat)

    plx_err = plxErr(G=G, vmini=V_IC, beta=ecl_lat)*1.0e-6

    return plx_err


def get_pos_err(G=None, V=None, V_IC=None, RA=None, Dec=None):
    """ A tool to determine the parallax uncertainty for a star

    Inputs
    ------
    G or V : float (provide one of these)
    Gaia G-band or Visual V-band magnitude of the star

    V_IC : float
    Color index of the stellar type

    RA : float
    Right ascension in degrees

    Dec: float
    Declination in degrees

    Returns
    -------
    pos_err : float
    position uncertainty in arcseconds

    """

    if V_IC is None:
        print("You must provide a color index in the form of V_IC.")
        return

    if G is None and V is None:
        print("You must provide at least one color: G or V")
        return

    if G is None:  G = V - 0.0257 - 0.0924 * (V_IC) - 0.1623 * (V_IC)**2 + 0.0090 * (V_IC)**3

    # KMB: get ecliptic latitude from ra/dec using astropy
    c = SkyCoord(ra=sys_ra*u.degree, dec=sys_dec*u.degree, frame='icrs')
    ecl_lat = c.GeocentricTrueEcliptic.lat

    # KMB: compute average N_obs for given ecliptic latitude
    N_obs_ave = N_transit_ave(ecl_lat)

    pos_err = posErr(G=G, vmini=V_IC, beta=ecl_lat)*1.0e-6

    return pos_err


def calc_new_magnitude(mag_current, distance_current, distance_new):
    """
    Calculate the magnitude of a star if its distance were shifted forward or backward

    Inputs
    ------
    mag_current : float
    Current magnitude of the stellar

    distance_current : float
    Current distance to the star (pc)

    distance_new : float
    New distance to the stellar (pc)

    Returns
    -------
    mag_new : float
    Magnitude of the star at the new distance

    """

    mag_new = mag_current + 5.0 * (np.log10(distance_new/10.0) - np.log10(distance_current/10.0))

    return mag_new
