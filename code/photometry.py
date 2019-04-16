import numpy as np
from scipy.interpolate import interp1d

# Global interpolation functions
M_G_interp = None
V_IC_interp = None

def load_stellar_photometry():
    """ Load stellar photometric data from file """

    global M_G_interp
    global V_IC_interp

    # Load up data
    data = np.genfromtxt("../data/mamajek_colors.dat", names=True)

    # Create V_IC interpolation object
    V_IC_interp = interp1d(data['Msun'], data['VIc'])

    # Create M_G interpolation object
    M_G = data['M_G']

    # M_G in data file does not go above 2 Msun - substitute with MV
    M_G[data['Msun']>2.0] = data['Mv'][data['Msun']>2.0]
    M_G_interp = interp1d(data['Msun'], M_G)

def get_M_G(mass):
    """
    Calculate the absolute G band magnitude

    Arguments
    ---------
    mass : float
        Input mass of a star

    Returns
    -------
    G_mag : float
        Return the G band magnitude
    """

    global M_G_interp

    if M_G_interp is None: load_stellar_photometry()
    return M_G_interp(mass)


def get_V_IC(mass):
    """
    Calculate the V_IC color

    Arguments
    ---------
    mass : float
        Input mass of a star

    Returns
    -------
    V_IC : float
        Return the V_IC color
    """

    global M_G_interp

    if M_G_interp is None: load_stellar_photometry()
    return M_G_interp(mass)


def get_G_mag_apparent(mass, distance):
    """
    Calculate the apparent G band magnitude

    Arguments
    ---------
    mass : float
        Stellar mass (Msun)
    distance : float
        Distance to the star (pc)

    Returns
    -------
    G_mag : float
        Apparent G band magnitude (mag)
    """
    return get_M_G(mass) + 5.0*np.log10(distance) - 5.0
