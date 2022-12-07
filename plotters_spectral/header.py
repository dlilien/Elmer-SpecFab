# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © N. M. Rathmann <rathmann@nbi.ku.dk>, 2022
#
# Distributed under terms of the GNU GPL3.0 license.


"""
Header file for common routines etc.
"""

import numpy as np
import scipy.special as sp

# This is the presumed psi_0^0 coefficient when deriving nlm from a^(2) or a^(4), a central normalization factor used below.
normfac = 1 / np.sqrt(4 * np.pi)
FS = 12


# fwd/inv transform pair for n_2^0 <--> a_zz^(2)
def n20_to_azz(x):
    return 1 / 3 * (1 + 2 / np.sqrt(5) * x)  # fwd


def azz_to_n20(x):
    return (x - 1 / 3) / (2 / np.sqrt(5))  # inv


def get_v_angles(a2, Ilami=-1):
    # Get a^(2) eigenvector angles, need to rotate (measured) ODFs into vertical frame needed for this analysis
    # If Ilami=-1, this will give the rotation angle for the largest eigenvalue direction (used for single max)
    # If Ilami=0, this will give the rotation angle for the smallest eigenvalue direction (used for girdles)
    lami, v = np.linalg.eig(a2)
    ind = np.argsort(lami)[Ilami]
    v1 = v[:, ind]
    v1_colat, v1_lon = np.arccos(v1[2]), np.arctan2(v1[1], v1[0])
    return (v1_colat, v1_lon, v1)


def get_deg(lat_or_colat, lon):
    return (np.rad2deg(lat_or_colat), np.rad2deg(lon))


def plot_trajectory(ax, nlm, arrpos=None, label=None, ls="-", lw=1.75, c="k", hwmul=1, zorder=14, endmarker=False, mse=7):

    x, y = np.real(nlm[:, 3]), np.real(nlm[:, 10])
    (h,) = ax.plot(x, y, ls=ls, lw=lw, c=c, label=label, zorder=zorder)

    if arrpos is not None:
        jj, djj = arrpos, 4
        xa0, ya0, xa1, ya1 = x[jj], y[jj], x[jj + djj], y[jj + djj]
        dx, dy = xa1 - xa0, ya1 - ya0
        arrlen = np.sqrt(dx**2 + dy**2)
        dx *= 0.5e-2 / arrlen
        dy *= 0.5e-2 / arrlen
        ax.arrow(xa0, ya0, dx, dy, shape="full", lw=1, color=c, head_width=0.02 * 2.2 * hwmul / normfac, head_length=0.035 * 0.8 / normfac, zorder=zorder)

    if endmarker:
        ax.plot(x[-1], y[-1], marker="o", ms=mse, ls="none", c=c, label=None, zorder=zorder)  # fillstyle='none',

    return h


def discretize_ODF(nlm, lm, latres=60):

    # latres = 60 # latitude resolution on S^2
    theta = np.linspace(0, np.pi, latres)  # CO-LAT
    phi = np.linspace(0, 2 * np.pi, 2 * latres)  # LON
    phi, theta = np.meshgrid(phi, theta)  # gridded
    lon, colat = phi, theta
    lat = np.pi / 2 - colat

    _, nlm_len = lm.shape
    F = np.real(np.sum([nlm[ii] * sp.sph_harm(lm[1][ii], lm[0][ii], phi, theta) for ii in np.arange(nlm_len)], axis=0))

    return (F, lon, lat)
