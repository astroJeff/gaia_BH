import numpy as np
from numpy.random import uniform
import emcee
from gatspy.periodic import LombScargleFast
from scipy.interpolate import interp1d
from pygaia.photometry import transformations
from pygaia.errors import spectroscopic

import prob
import orbit
import gaia_tools as gaia
import const as c
from astropy import units as u
from astropy.coordinates import SkyCoord
import plotting


def load_data(filename):
    dtype = [("M1","<f8"), ("M2","<f8"), ("ecc","<f8"), ("P_orb","<f8"),
             ("Lum2", "<f8"), ("Temp2", "<f8"), ("Rad2", "<f8"),
             ("Xgx", "<f8"), ("Ygx", "<f8"), ("Zgx", "<f8"), ("dist","<f8"),
             ("inc","<f8"), ("Omega","<f8"), ("omega","<f8"),
             ("proj_sep","<f8"), ("G","<f8"), ("V_IC","<f8")]

    data = np.genfromtxt(filename, dtype=dtype, delimiter=',')

    return data



def get_M2_RV_err(mass, G, V_IC):

    # First, get a V magnitude
    G_VC = transformations.gminvFromVmini(V_IC)
    V = G - G_VC

    # Adjust for masses below pyGaia limit
    if mass < 0.67: return spectroscopic.vradErrorSkyAvg(V, 'K4V')
    if mass > 17.5: return spectroscopic.vradErrorSkyAvg(V, 'B0V')

    # Stellar masses from Carol & Ostlie
    # To Do: Check these numbers
    mass_grid = np.array([17.5, 5.9, 2.9, 2.0, 1.6, 1.05, 1.0, 0.79, 0.67])
    str_grid = np.array(['B0V','B5V','A0V','A5V','F0V','G0V','G5V','K0V','K4V'])
    rv_err_grid = np.array([])
    for i, m in enumerate(mass_grid):
        rv_err_tmp = spectroscopic.vradErrorSkyAvg(V, str_grid[i])
        rv_err_grid = np.append(rv_err_grid, rv_err_tmp)

    rv_err_interp = interp1d(mass_grid, rv_err_grid)

    return rv_err_interp(mass)



def get_truths(sys):

    coords = SkyCoord(x=sys['Xgx']*u.parsec, y=sys['Ygx']*u.parsec, z=sys['Zgx']*u.parsec,
                      frame='galactocentric')

    sys_ra = coords.icrs.ra.deg
    sys_dec = coords.icrs.dec.deg

    I = sys['inc']
    Omega = sys['Omega']
    omega = sys['omega']
    e = sys['ecc']
    P = sys['P_orb'] * c.secyer
    M1 = sys['M1'] * c.Msun
    M2 = sys['M2'] * c.Msun
    distance = sys['dist'] * 1.0e3
    tau = 1.0e4
    gamma = 10.0

    truths = sys_ra, sys_dec, Omega, omega, I, tau, e, P, gamma, M1, M2, distance

    return truths


def run_one_binary(i, sys, nwalkers=100, ntemps=1, nburn=1000, nsteps=2000):

    # create array from input system parameters
    obs_pos = create_data_array(sys)

    # Set companion mass constraint
    M2_obs_err = 0.5
    M2_obs = np.random.normal(loc=sys['M2'], scale=M2_obs_err, size=1)

    # Calculate truths for the binary
    truths = get_truths(sys)

    # Initialize the emcee sampler
    p0 = initialize_sampler(truths, obs_pos, nwalkers=nwalkers, ntemps=ntemps)

    # Run the sampler
    chains = run_emcee(p0, obs_pos, M2_obs, M2_obs_err, nburn=nburn, nsteps=nsteps)

    # Save the chains
    file_out = "../data/binary_set/system_" + "%.i" % i + "_chains.npy"
    np.save(file_out, chains)





    # Plot data
    if chains.ndim == 4:  # For PT sampler
        chains_in = chains[0]
    else:
        chains_in = chains

    nchains, nsteps, nvar = chains_in.shape
    chains_flat = chains_in.reshape((nchains*nsteps, nvar))

    plotting.plot_corner(chains_flat, truths, filename="../figures/binary_set/binary_" + "%.i" % i + "_corner.pdf")
    plotting.plot_masses(chains_flat, truths, filename="../figures/binary_set/binary_" + "%.i" % i + "_masses.pdf")
    plotting.plot_orbit(obs_pos, chains_flat, truths, filename="../figures/binary_set/binary_" + "%.i" % i + "_orbit.pdf")
    plotting.plot_trace(chains_in, filename="../figures/binary_set/binary_" + "%.i" % i + "_trace.pdf")


def create_data_array(sys):

    names = ["time","ra","dec","ra_err","dec_err","rv","rv_err","plx","plx_err"]
    obs_pos = np.recarray(c.N_samples, names=names,
                          formats=['float64,float64,float64,float64,float64,float64,float64,float64,float64'])

    time = np.sort(uniform(c.time_max*c.secyer, size=c.N_samples))


    # save system parameters
    p = get_truths(sys)


    # Calculate orbit
    ra, dec = orbit.get_ra_dec(p, time)
    rv = orbit.get_RV(p, time)

    # Save data to array
    obs_pos['time'] = time  # in seconds
    obs_pos['ra'] = ra  # in degrees
    obs_pos['dec'] = dec  # in degrees
    obs_pos['rv'] = rv  # in km/s

    obs_pos['ra_err'] = gaia.get_single_obs_pos_err(G=sys[15], V_IC=sys[16], RA=ra, Dec=dec,
                                          DIST=sys[10], XGX=sys[7], YGX=sys[8], ZGX=sys[9]) # in arcseconds
    obs_pos['dec_err'] = gaia.get_single_obs_pos_err(G=sys[15], V_IC=sys[16], RA=ra, Dec=dec,
                                          DIST=sys[10], XGX=sys[7], YGX=sys[8], ZGX=sys[9]) # in arcseconds
    # obs_pos['ra_err'] = 5.34e-6 / 3600.0 # in deg
    # obs_pos['dec_err'] = 5.34e-6 / 3600.0 # in deg

    obs_pos['rv_err'] = get_M2_RV_err(sys['M2'], sys['G'], sys['V_IC'])  # in km/s

    # parallax
    obs_pos['plx'] = 1.0/(sys['dist']*1.0e3) # plx in arseconds
    obs_pos['plx_err'] = gaia.get_plx_err(G=sys[15], V_IC=sys[16], RA=ra, Dec=dec,
                                          DIST=sys[10], XGX=sys[7], YGX=sys[8], ZGX=sys[9]) # in arcseconds
    # obs_pos['plx_err'] = gaia.get_plx_err(G=sys[10], V_IC=0.75)*1.0e3  # We use the V_IC of a G2V star

    # Add uncertainties to measurements
    obs_pos['ra'] += np.random.normal(scale=obs_pos['ra_err'], size=len(obs_pos))
    obs_pos['dec'] += np.random.normal(scale=obs_pos['dec_err'], size=len(obs_pos))
    obs_pos['rv'] += np.random.normal(scale=obs_pos['rv_err'], size=len(obs_pos))

    return obs_pos


def initialize_sampler(truths, obs_pos, nwalkers=100, ntemps=1):

    ndim = 12

    # Starting array
    if ntemps == 1:
        p0 = np.zeros((nwalkers, ndim))
        delta = 1.0 + np.random.normal(0.0, 0.001, size=(nwalkers, ndim)) # Random deviations from best fit
        for i in range(ndim): p0[:,i] = truths[i] * delta[:,i] # Adjust starting array accordingly

    else:
        p0 = np.zeros((ntemps, nwalkers, ndim))
        delta = 1.0 + np.random.normal(0.0, 0.001, size=(ntemps, nwalkers, ndim)) # Random deviations from best fit
        for i in range(ndim): p0[:,:,i] = truths[i] * delta[:,:,i] # Adjust starting array accordingly


    return p0


def run_emcee(p0, obs_pos, M2_obs, M2_obs_err, nburn=1000, nsteps=2000, thin=10):

    if p0.ndim == 2:
        nwalkers, ndim = p0.shape
        ntemps = 1
    elif p0.ndim == 3:
        ntemps, nwalkers, ndim = p0.shape


    # Define sampler
    if ntemps == 1:
        sampler = emcee.EnsembleSampler(nwalkers=nwalkers, dim=ndim, lnpostfn=prob.ln_posterior, args=(obs_pos, M2_obs, M2_obs_err))
    else:
        sampler = emcee.PTSampler(ntemps=ntemps, nwalkers=nwalkers, dim=ndim,
                                  logl=prob.ln_likelihood, logp=prob.ln_prior,
                                  loglargs=(obs_pos, M2_obs, M2_obs_err))


    # Run emcee
    if ntemps == 1:
        # Run burn-in
        pos,ln_prob,state = sampler.run_mcmc(p0, N=nburn)

        # Run sampler
        sampler.reset()
        pos,ln_prob,state = sampler.run_mcmc(pos, N=nsteps)

        # Downsample by factor of 10
        chains = sampler.chain[:,::thin,:]

    else:  # Parallel tempering algorithm
        # Run burn-in
        for pos,ln_prob,state in sampler.sample(p0, iterations=nburn):  pass

        # Run sampler
        sampler.reset()
        for pos,ln_prob,state in sampler.sample(pos, iterations=nsteps, thin=thin):  pass

        # return chains
        chains = sampler.chain

    return chains


def run_set():

    filename = "../data/GxExampleSorted_m2_short.dat"

    data = load_data(filename)

    print("Cycling through binaries...")

    for i, sys in enumerate(data):

        if i != 0: continue

        print("...running", i, "of", len(data), "binaries")

        run_one_binary(i, sys, nwalkers=100, ntemps=1, nburn=200, nsteps=500)


run_set()
