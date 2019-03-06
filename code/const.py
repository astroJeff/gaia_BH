import numpy as np

# Constants in cgs
Msun = 1.989e33
G = 6.674e-8
secday = 3600.0*24.0
secyer = 3600.0*24.0*365.25
AUincm = 1.496e13
AUinm = 1.4959787066e11
pcincm = 3.086e18

rad_to_deg = 180.0 / np.pi
deg_to_rad = np.pi / 180.0
J2000 = 2451545.0  # Julian date definition
FLATTEN = 0.003352813
EQUAT_RAD = 6378137.   # equatorial radius of earth, meters

# Gaia numbers
time_max = 5.0  # in years
N_samples = 70  # over 5 year lifetime
