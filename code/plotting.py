import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import corner

import const as c
import orbit




def plot_pos(p, obs_pos):

    t_tmp = np.linspace(np.min(obs_pos['time']), np.max(obs_pos['time']), 1000)

    if type(p) is tuple or p.ndim==1:
        if len(p)==8:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], 0.0, p[4], p[5], p[6], p[7], 1.0e3/obs_pos['plx'][0]
        elif len(p)==9:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], 1.0e3/obs_pos['plx'][0]
        else:
            p_full = p

        ra, dec = orbit.get_ra_dec(p_full, t_tmp)
        plt.plot(ra*3600.0*1.0e3, dec*3600.0*1.0e3, color='k')
    else:
        for x in p:
            if len(x)==8:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], 0.0, x[4], x[5], x[6], x[7], 1.0e3/obs_pos['plx'][0]
            elif len(x)==9:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], 1.0e3/obs_pos['plx'][0]
            else:
                p_full = x

            ra, dec = orbit.get_ra_dec(p_full, t_tmp)
            plt.plot(ra*3600.0*1.0e3, dec*3600.0*1.0e3, alpha=0.05, color='k')


    plt.scatter(obs_pos['ra']*3600.0*1.0e3, obs_pos['dec']*3600.0*1.0e3)

#     plt.xlim(-2, 2)
#     plt.ylim(-2, 2)
    plt.xlabel(r"$\Delta\ \alpha$ (mas)")
    plt.ylabel(r"$\Delta\ \delta$ (mas)")
    plt.axis('equal')

    plt.show()

def plot_pos_ra(p, obs_pos):

    t_tmp = np.linspace(np.min(obs_pos['time']), np.max(obs_pos['time']), 1000)

    if type(p) is tuple or p.ndim==1:
        if len(p)==8:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], 0.0, p[4], p[5], p[6], p[7], 1.0e3/obs_pos['plx'][0]
        elif len(p)==9:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], 1.0e3/obs_pos['plx'][0]
        else:
            p_full = p

        ra, dec = orbit.get_ra_dec(p_full, t_tmp)
        plt.plot(t_tmp/c.secyer, ra*3600.0*1.0e3, color='k')
    else:
        for x in p:
            if len(x)==8:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], 0.0, x[4], x[5], x[6], x[7], 1.0e3/obs_pos['plx'][0]
            elif len(x)==9:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], 1.0e3/obs_pos['plx'][0]
            else:
                p_full = x

            ra, dec = orbit.get_ra_dec(p_full, t_tmp)
            plt.plot(t_tmp/c.secyer, ra*3600.0*1.0e3, alpha=0.05, color='k')


    ra_tmp, dec_tmp = orbit.get_ra_dec(p_full, obs_pos['time'])

    plt.scatter(obs_pos['time']/c.secyer, obs_pos['ra']*3600.0*1.0e3)

#     plt.ylim(-2, 2)
    plt.xlabel("Time (yr)")
    plt.ylabel(r"$\Delta\ \alpha$ (mas)")

    plt.show()

def plot_pos_dec(p, obs_pos):

    t_tmp = np.linspace(np.min(obs_pos['time']), np.max(obs_pos['time']), 1000)


    if type(p) is tuple or p.ndim==1:
        if len(p)==8:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], 0.0, p[4], p[5], p[6], p[7], 1.0e3/obs_pos['plx'][0]
        elif len(p)==9:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], 1.0e3/obs_pos['plx'][0]
        else:
            p_full = p

        ra, dec = orbit.get_ra_dec(p_full, t_tmp)
        plt.plot(t_tmp/c.secyer, dec*3600.0*1.0e3, color='k')
    else:
        for x in p:
            if len(x)==8:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], 0.0, x[4], x[5], x[6], x[7], 1.0e3/obs_pos['plx'][0]
            elif len(x)==9:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], 1.0e3/obs_pos['plx'][0]
            else:
                p_full = x

            ra, dec = orbit.get_ra_dec(p_full, t_tmp)
            plt.plot(t_tmp/c.secyer, dec*3600.0*1.0e3, alpha=0.05, color='k')


    plt.scatter(obs_pos['time']/c.secyer, obs_pos['dec']*3600.0*1.0e3)

#     plt.ylim(-2, 2)
    plt.xlabel("Time (yr)")
    plt.ylabel(r"$\Delta\ \delta$ (mas)")

    plt.show()


def plot_rv(p, obs_pos):

    t_tmp = np.linspace(np.min(obs_pos['time']), np.max(obs_pos['time']), 1000)


    if type(p) is tuple or p.ndim==1:
        if len(p)==8:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], 0.0, p[4], p[5], p[6], p[7], 1.0e3/obs_pos['plx'][0]
        elif len(p)==9:
            p_full = 0.0, 0.0, p[0], p[1], p[2], p[3], p[4], p[5], p[6], p[7], p[8], 1.0e3/obs_pos['plx'][0]
        else:
            p_full = p

        rv = orbit.get_RV(p_full, t_tmp)
        plt.plot(t_tmp/c.secyer, rv, color='k')
    else:
        for x in p:
            if len(x)==8:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], 0.0, x[4], x[5], x[6], x[7], 1.0e3/obs_pos['plx'][0]
            elif len(x)==9:
                p_full = 0.0, 0.0, x[0], x[1], x[2], x[3], x[4], x[5], x[6], x[7], x[8], 1.0e3/obs_pos['plx'][0]
            else:
                p_full = x

            rv = orbit.get_RV(p_full, t_tmp)
            plt.plot(t_tmp/c.secyer, rv, alpha=0.05, color='k')


    plt.scatter(obs_pos['time']/c.secyer, obs_pos['rv'])
    plt.xlabel("Time (yr)")
    plt.ylabel("Radial velocity (km s$^{-1}$)")

    plt.show()


def plot_trace(chains, alpha=0.05, filename=None):

    fig, ax = plt.subplots(6, 2, figsize=(12, 18))

    labels=[r'$\alpha$', r'$\delta$', r'$\Omega$', r'$\omega$', r'$i$', r'$\tau$', r'$e$', r'$P_{\rm orb}$',
            r'$\gamma$', r'M$_1$', r'M$_2$', 'distance']

    nwalker, nstep, ndim = chains.shape

    for j in range(ndim):

        jx = int(j/2)
        jy = j%2

        for i in range(nwalker):
            ax[jx,jy].plot(chains[i,:,j], color='k', alpha=alpha, rasterized=True)

        ax[jx,jy].set_ylabel(labels[j])

        ax[jx,jy].set_xlabel("Step")


    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename, rasterized=True)
    else:
        plt.show()

    plt.clf()


def plot_corner(chains, truths, filename=None):
    labels=[r'$\alpha$', r'$\delta$', r'$\Omega$', r'$\omega$', r'$i$', r'$\tau$', r'$e$', r'$P_{\rm orb}$', r'$\gamma$', \
            r'M$_1$', r'M$_2$', 'distance']

    corner.corner(chains, truths=truths, labels=labels)

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

    plt.clf()


def plot_masses(chains, truths, filename=None):

    plt.figure(figsize=(6,6))

    names = [r"M$_1$ (M$_{\odot}$)", r"M$_2$ (M$_{\odot}$)"]

    size = len(chains[:,9])
    M1_sorted = np.sort(chains[:,9])
    M2_sorted = np.sort(chains[:,10])


    lims = [(np.min(M1_sorted/c.Msun), M1_sorted[int(0.99*size)]/c.Msun),
            (np.min(M2_sorted/c.Msun), M2_sorted[int(0.99*size)]/c.Msun)]
#     lims = [(0.9*np.min(chains[:,9]/c.Msun), 1.1*np.max(chains[:,9]/c.Msun)),
#             (0.9*np.min(chains[:,10]/c.Msun), 1.1*np.max(chains[:,10]/c.Msun))]
#     lims=[(6,14), (16,24)]
#     print(lims, truths[9]/c.Msun, truths[10]/c.Msun)
    inputs = [truths[9]/c.Msun, truths[10]/c.Msun]

#     plt.scatter(chains[:,9]/c.Msun, chains[:,10]/c.Msun, alpha=0.1)

    corner.hist2d(chains[:,9]/c.Msun, chains[:,10]/c.Msun, range=lims, bins=40)

    plt.xlabel(names[0])
    plt.ylabel(names[1])
    plt.axhline(truths[10]/c.Msun, color='C0')
    plt.axvline(truths[9]/c.Msun, color='C0')

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

    plt.clf()


def plot_orbit(obs_pos, chains, truths, filename=None, xlim=None,
               nsamples=30, alpha=0.1):

    nchains, nsteps, nvar = chains.shape
    chains_flat = chains.reshape((nchains*nsteps, nvar))
    idx = np.random.randint(nchains*nsteps, size=nsamples)


    fig, ax = plt.subplots(2, 2, figsize=(8,8))

    ax[0,0].scatter(obs_pos['ra']*3600.0*1.0e3, obs_pos['dec']*3600.0*1.0e3)

    t_tmp = np.linspace(0.0, np.max(obs_pos['time']), 10000)
    ra, dec = orbit.get_ra_dec(truths, t_tmp)
    # ax[0,0].plot(ra*3600.0*1.0e3, dec*3600.0*1.0e3, color='k')
    for i in range(nsamples):
        ra, dec = orbit.get_ra_dec(chains_flat[idx[i]], t_tmp)
        ax[0,0].plot(ra*3600.0*1.0e3, dec*3600.0*1.0e3, color='k', alpha=alpha, rasterized=True)


    ax[0,0].set_xlabel(r"$\Delta\ \alpha$ (mas)")
    ax[0,0].set_ylabel(r"$\Delta\ \delta$ (mas)")


    rv = orbit.get_RV(truths, t_tmp)
    # ax[1,0].plot(t_tmp/c.secyer, rv, color='k')
    for i in range(nsamples):
        rv = orbit.get_RV(chains_flat[idx[i]], t_tmp)
        ax[1,0].plot(t_tmp/c.secyer, rv, color='k', alpha=alpha)

    ax[1,0].scatter(obs_pos['time']/c.secyer, obs_pos['rv'], zorder=100)
    ax[1,0].set_xlabel("Time (yr)")
    ax[1,0].set_ylabel("Radial velocity (km s$^{-1}$)")
    if xlim is not None: ax[1,0].set_xlim(xlim)


    # ax[0,1].plot(t_tmp/c.secyer, ra*3600.0*1.0e3, color='k')
    for i in range(nsamples):
        ra, dec = orbit.get_ra_dec(chains_flat[idx[i]], t_tmp)
        ax[0,1].plot(t_tmp/c.secyer, ra*3600.0*1.0e3, color='k', alpha=alpha)

    ax[0,1].scatter(obs_pos['time']/c.secyer, obs_pos['ra']*3600.0*1.0e3, zorder=100)
    ax[0,1].set_xlabel("Time (yr)")
    ax[0,1].set_ylabel(r"$\Delta\ \alpha$ (mas)")

    P_orb = truths[7]/c.secyer
    ax[0,1].set_title(r"P$_{\rm orb}$ = " + "%.4f" % P_orb + " yr")
    if xlim is not None: ax[0,1].set_xlim(xlim)


    # ax[1,1].plot(t_tmp/c.secyer, dec*3600.0*1.0e3, color='k')
    for i in range(nsamples):
        ra, dec = orbit.get_ra_dec(chains_flat[idx[i]], t_tmp)
        ax[1,1].plot(t_tmp/c.secyer, dec*3600.0*1.0e3, color='k', alpha=alpha)

    ax[1,1].scatter(obs_pos['time']/c.secyer, obs_pos['dec']*3600.0*1.0e3, zorder=100)
    ax[1,1].set_xlabel("Time (yr)")
    ax[1,1].set_ylabel(r"$\Delta\ \delta$ (mas)")
    if xlim is not None: ax[1,1].set_xlim(xlim)

    plt.tight_layout()

    if filename is not None:
        plt.savefig(filename)
    else:
        plt.show()

    plt.clf()
