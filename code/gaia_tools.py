import numpy as np


def get_single_obs_err(G=None, V=None, V_IC=None):
    """ A tool to determine the single observation astrometric precision. This is
    calculated in a bit of an ad hoc way. We call the overall parallax precision, then
    multiply by sqrt(70), as this is the number that should have been used to determine
    the parallax precision.

    Inputs
    ------
    G or V : float (provide one of these)
    Gaia G-band or Visual V-band magnitude of the star

    V_IC : float
    Color index of the stellar type

    Returns
    -------
    sigma_plx : float
    parallax uncertainty in arcseconds

    """


    if V_IC is None:
        print("You must provide a color index in the form of V_IC.")
        return

    if G is None and V is None:
        print("You must provide at least one color: G or V")
        return


    plx_err = get_plx_err(G=G, V=V, V_IC=V_IC)

    m = 1.2
    g_pi = 1.9

    pos_err = np.sqrt(70.0) / (m * g_pi) * plx_err

    return pos_err




def get_plx_err(G=None, V=None, V_IC=None):
    """ A tool to determine the parallax uncertainty for a star

    Inputs
    ------
    G or V : float (provide one of these)
    Gaia G-band or Visual V-band magnitude of the star

    V_IC : float
    Color index of the stellar type

    Returns
    -------
    sigma_plx : float
    parallax uncertainty in arcseconds

    """


    if V_IC is None:
        print("You must provide a color index in the form of V_IC.")
        return

    if G is None and V is None:
        print("You must provide at least one color: G or V")
        return

    if G is None:  G = V - 0.0257 - 0.0924 * (V_IC) - 0.1623 * (V_IC)**2 + 0.0090 * (V_IC)**3

    z = np.max([10**(0.4 * (12.09 - 15.0)), 10**(0.4 * (G - 15.0))])
    sigma_plx = np.sqrt(-1.631 + 680.766*z + 32.732*z**2) * (0.986 + (1 - 0.986) * V_IC)

    # Convert from microas to asec
    sigma_plx = sigma_plx/1.0e6

    return sigma_plx



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
