import numpy as np
from numpy.random import uniform
import emcee
from gatspy.periodic import LombScargleFast

import prob
import orbit
import gaia_tools as gaia
import const as c
import plotting


def load_data(filename):
    dtype = [("M1","<f8"), ("M2","<f8"), ("ecc","<f8"), ("P_orb","<f8"), ("sep","<f8"),
         ("proj_sep","<f8"), ("dist","<f8"), ("inc","<f8"), ("Omega","<f8"),
         ("omega","<f8"), ("G","<f8"), ("sigma","<f8")]

    data = np.genfromtxt(filename, dtype=dtype, delimiter=',')

    return data


def get_truths(sys):

    sys_ra = 100.0
    sys_dec = 10.0
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

    # Calculate truths for the binary
    truths = get_truths(sys)

    # Initialize the emcee sampler
    p0 = initialize_sampler(truths, obs_pos, nwalkers=nwalkers, ntemps=ntemps)

    # Run the sampler
    chains = run_emcee(p0, obs_pos, nburn=nburn, nsteps=nsteps)

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

    obs_pos['ra_err'] = 1.0e-4 / 3600.0 # in deg
    obs_pos['dec_err'] = 1.0e-4 / 3600.0 # in deg
    # obs_pos['ra_err'] = 5.34e-6 / 3600.0 # in deg
    # obs_pos['dec_err'] = 5.34e-6 / 3600.0 # in deg
    obs_pos['rv_err'] = 1.0  # in km/s

    # parallax
    obs_pos['plx'] = 1.0/(sys['dist']*1.0e3) # plx in arseconds
    obs_pos['plx_err'] = sys['sigma']  # Use the parallax errors provided by Katie
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


def run_emcee(p0, obs_pos, nburn=1000, nsteps=2000, thin=10):

    if p0.ndim == 2:
        nwalkers, ndim = p0.shape
        ntemps = 1
    elif p0.ndim == 3:
        ntemps, nwalkers, ndim = p0.shape


    # Define sampler
    if ntemps == 1:
        sampler = emcee.EnsembleSampler(nwalkers=nwalkers, dim=ndim, lnpostfn=prob.ln_posterior, args=(obs_pos,))
    else:
        sampler = emcee.PTSampler(ntemps=ntemps, nwalkers=nwalkers, dim=ndim,
                                  logl=prob.ln_likelihood, logp=prob.ln_prior,
                                  loglargs=(obs_pos,))


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

    filename = "../data/gxSortedByDistance.dat"

    data = load_data(filename)

    print("Cycling through binaries...")

    for i, sys in enumerate(data):

        if i != 2: continue

        print("...running", i, "of", len(data), "binaries")

        run_one_binary(i, sys, nwalkers=100, ntemps=1, nburn=2000, nsteps=5000)


# run_set()
