import numpy as np

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

import os

Msun = 1.989e33

data_dir = '../data/binary_chains/'

dist_list = np.array([])
rv_err_list = np.array([])


# Determine grid size
for file_name in os.listdir(data_dir):

    elements = file_name.split('_')

    distance = elements[0]
    distance = float(distance[8:])

    rv_err = elements[1]
    rv_err = rv_err[5:]
    rv_err = float(rv_err[:-4])

    if not distance in dist_list: dist_list = np.append(dist_list, distance)
    if not rv_err in rv_err_list: rv_err_list = np.append(rv_err_list, rv_err)


# Sort the distance and rv_err set
dist_list = np.sort(dist_list)
rv_err_list = np.sort(rv_err_list)


# Get the lengths of the two axes
n_dist = len(dist_list)
n_rv_err = len(rv_err_list)


# Mass errors
M1_err = np.zeros((n_dist, n_rv_err))
M2_err = np.zeros((n_dist, n_rv_err))


# Load data
for file_name in os.listdir(data_dir):

    elements = file_name.split('_')

    distance = elements[0]
    distance = float(distance[8:])

    rv_err = elements[1]
    rv_err = rv_err[5:]
    rv_err = float(rv_err[:-4])

    data = np.load(os.path.join(data_dir, file_name))

    M1 = data[:,:,9]
    M2 = data[:,:,10]

    n_walkers, n_steps = M1.shape

    M1 = M1.reshape(n_walkers*n_steps)
    M2 = M2.reshape(n_walkers*n_steps)

    M1_sorted = np.sort(M1)
    M2_sorted = np.sort(M2)

    idx1 = int(0.13*n_walkers*n_steps)
    idx2 = int(0.87*n_walkers*n_steps)

    idx_dist = np.where(distance == dist_list)[0]
    idx_rv_err = np.where(rv_err == rv_err_list)[0]

    M1_err[idx_dist, idx_rv_err] = 0.5*(M1_sorted[idx2] - M1_sorted[idx1])/Msun
    M2_err[idx_dist, idx_rv_err] = 0.5*(M2_sorted[idx2] - M2_sorted[idx1])/Msun



# Plot grid
X, Y = np.meshgrid(dist_list, rv_err_list)

im = plt.tricontourf(X.flatten(), Y.flatten(), M1_err.T.flatten(), cmap='Blues')
plt.colorbar(im)

plt.xscale('log')
plt.xlim(50, 10000)

plt.yscale('log')
plt.ylim(0.3, 20.0)

plt.xlabel("Distance (pc)")
plt.ylabel("RV Precision (km s$^{-1}$)")

plt.title("Mass Uncertainty (M$_{\odot}$)")

plt.tight_layout()
plt.savefig("../figures/M1_grid.pdf")

    # print(distance, rv_err, M1_err/Msun, M2_err/Msun)
