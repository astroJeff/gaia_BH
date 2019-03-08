import numpy as np
from scipy.interpolate import interp1d


M_G_interp = None
V_IC_interp = None

def load_stellar_photometry():

    global M_G_interp
    global V_IC_interp

    data = np.genfromtxt("../data/mamajek_colors.dat", names=True)
    V_IC_interp = interp1d(data['Msun'], data['VIc'])
    M_G = data['M_G']
    M_G[data['Msun']>2.0] = data['Mv'][data['Msun']>2.0]
    M_G_interp = interp1d(data['Msun'], M_G)

def get_M_G(mass):

    global M_G_interp

    if M_G_interp is None: load_stellar_photometry()
    return M_G_interp(mass)


def get_V_IC(mass):

    global M_G_interp

    if M_G_interp is None: load_stellar_photometry()
    return M_G_interp(mass)

def get_G_mag_apparent(mass, distance):
    return get_M_G(mass) + 5.0*np.log10(distance) - 5.0
