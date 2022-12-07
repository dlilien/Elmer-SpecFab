# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © D. A. Lilien <david.lilien@umanitoba.ca>, 2022
#
# Distributed under terms of the GNU GPL3.0 license.


import numpy as np
import scipy.special as sp
import matplotlib.pyplot as plt
from matplotlib.markers import MarkerStyle
import matplotlib.gridspec as gridspec
import cartopy.crs as ccrs

# from .anisolib import cartoon_norm

from specfabpy import specfabpy as sf

inclination = 60
rot = -90 * 1
PRJ = ccrs.Orthographic(rot, 90 - inclination)

HOOP = [(0.25 * np.cos(x), np.sin(x)) for x in np.linspace(0, 2 * np.pi, 50)]


def dumb_woodcock(stuff_dict):
    norm = np.abs(stuff_dict["eigenv 3"]) + np.abs(stuff_dict["eigenv 2"]) + np.abs(stuff_dict["eigenv 1"])
    e1 = np.abs(stuff_dict["eigenv 1"]) / norm
    e2 = np.abs(stuff_dict["eigenv 2"]) / norm
    e3 = np.abs(stuff_dict["eigenv 3"]) / norm
    out = e3 / e2
    out[e1 == 0] = 1e8
    out[e2 == 0] = 1e8
    return out > 5.0


def woodcock(stuff_dict):
    norm = np.abs(stuff_dict["eigenv 3"]) + np.abs(stuff_dict["eigenv 2"]) + np.abs(stuff_dict["eigenv 1"])
    e1 = np.abs(stuff_dict["eigenv 1"]) / norm
    e2 = np.abs(stuff_dict["eigenv 2"]) / norm
    e3 = np.abs(stuff_dict["eigenv 3"]) / norm
    out = np.log(e3 / e2) / np.log(e2 / e1)
    return out


def a2plot(a2, ax=None, cbr=True, show=True, latres=20, lonres=40, levels=8, Omax=0.3, Omin=0.0, cmap="Greys"):
    phi = np.linspace(0, 2 * np.pi, lonres)  # LON
    theta = np.linspace(0, np.pi, latres)  # CO-LAT
    phi, theta = np.meshgrid(phi, theta)
    lon, colat = phi, theta
    lat = np.pi / 2 - colat

    x = 0
    y = 1
    z = 2
    k = np.sqrt(15.0 / (2.0 * np.pi))

    # ODF expansion coefs from a^2
    c00 = (a2[x, x] + a2[y, y] + a2[z, z]) / np.sqrt(4 * np.pi)
    c20 = -1 / 4 * np.sqrt(5 / np.pi) * (a2[x, x] + a2[y, y] - 2 * a2[z, z])
    c22m = k / 4 * (a2[x, x] + 2j * a2[x, y] - a2[y, y])  # c_2^{-2}
    c22p = k / 4 * (a2[x, x] - 2j * a2[x, y] - a2[y, y])  # c_2^{+2}
    c21m = +k / 2 * (a2[x, z] + 1j * a2[y, z])
    c21p = -k / 2 * (a2[x, z] - 1j * a2[y, z])

    # Discretize ODF
    lm = [(0, 0), (2, -2), (2, -1), (2, 0), (2, 1), (2, 2)]
    clm = [c00, c22m, c21m, c20, c21p, c22p]
    ODF = np.sum([clm[ii] * sp.sph_harm(m, l, phi, theta) for ii, (l, m) in enumerate(lm)], axis=0)
    ODF = np.real(ODF)  # Any imag values are rounding errors.
    ODF[ODF < 0.0] = 0.0

    if ax is None:
        plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax = plt.subplot(gs[0, 0], projection=PRJ)
        ax.set_global()
    # cmap.set_under('#fee0d2')
    # ODF[ODF < 0.0] = 0.0001
    # h1 = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), ODF, transform=ccrs.PlateCarree(), norm=mpl.colors.LogNorm(0.01, 1), levels=np.logspace(-2, 0, 11), extend='min', cmap="Greys")
    if Omin < 0:
        extend = "both"
    else:
        extend = "max"
    h1 = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), ODF, transform=ccrs.PlateCarree(), levels=np.linspace(Omin, Omax, levels), extend=extend, cmap=cmap)
    ax.gridlines(crs=ccrs.PlateCarree())
    if cbr:
        cb1 = plt.colorbar(h1, ax=ax, fraction=0.035, aspect=19, orientation="horizontal", pad=0.1)
        cb1.set_label(r"ODF$(\theta,\phi)$ - mean(ODF) (i.e. isotropic contr. removed)")
    if show:
        print(np.mean(ODF))
        plt.show()
    return h1


def a4plot(a2, a4, elmer=True, ax=None, cbr=True, show=True, latres=20, lonres=40, levels=8, Omax=0.3):
    """Modified from Nicholas' specfab version--add parameters so it can go in other plots"""
    phi = np.linspace(0, 2 * np.pi, lonres)  # LON
    theta = np.linspace(0, np.pi, latres)  # CO-LAT
    phi, theta = np.meshgrid(phi, theta)
    lon, colat = phi, theta
    lat = np.pi / 2 - colat

    # Determine spectral coefs from a^(2) and a^(4)
    lm, nlm_len = sf.init(4)

    if elmer:
        clm = sf.get_nlm_elmer_a4(nlm_len, a2, a4)
    else:
        clm = sf.get_a4_to_nlm(nlm_len, a2, a4)

    ODF = np.sum([clm[ii] * sp.sph_harm(m, l, phi, theta) for ii, (l, m) in enumerate(lm.T)], axis=0)
    ODF = np.real(ODF)  # Any imag values are rounding errors.

    if ax is None:
        plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax = plt.subplot(gs[0, 0], projection=PRJ)
        ax.set_global()
    ODF[ODF < 0.0] = 0.0
    h1 = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), ODF, transform=ccrs.PlateCarree(), levels=np.linspace(0, Omax, levels), extend="max", cmap="Greys")
    ax.gridlines(crs=ccrs.PlateCarree())
    if cbr:
        cb1 = plt.colorbar(h1, ax=ax, fraction=0.035, aspect=19, orientation="horizontal", pad=0.1)
        cb1.set_label(r"ODF$(\theta,\phi)$ - mean(ODF) (i.e. isotropic contr. removed)")
    if show:
        print(np.mean(ODF))
        plt.show()
    return h1


def nlm_plot(LCap, c, elmer=True, ax=None, cbr=True, show=True, latres=20, lonres=40, levels=8, Omin=0.0, Omax=0.3, cmap="Greys"):
    """Modified from Nicholas' specfab version--add parameters so it can go in other plots"""
    phi = np.linspace(0, 2 * np.pi, lonres)  # LON
    theta = np.linspace(0, np.pi, latres)  # CO-LAT
    phi, theta = np.meshgrid(phi, theta)
    lon, colat = phi, theta
    lat = np.pi / 2 - colat

    # Determine spectral coefs from a^(2) and a^(4)
    lm, nlm_len = sf.init(LCap)

    if elmer:
        clm = c[:nlm_len] + 1.0j * c[nlm_len:]
    else:
        clm = c

    ODF = np.sum([clm[ii] * sp.sph_harm(m, l, phi, theta) for ii, (l, m) in enumerate(lm.T)], axis=0)
    ODF = np.real(ODF)  # Any imag values are rounding errors.
    ODF = np.divide(ODF, clm[0] * np.sqrt(4.0 * np.pi))
    ODF[ODF < 0.0] = 0.0

    if ax is None:
        plt.figure()
        gs = gridspec.GridSpec(1, 1)
        ax = plt.subplot(gs[0, 0], projection=PRJ)
        ax.set_global()
    # ODF[ODF < 0.0] = 0.0
    if Omin < 0:
        extend = "both"
    else:
        extend = "max"
    h1 = ax.contourf(np.rad2deg(lon), np.rad2deg(lat), ODF, transform=ccrs.PlateCarree(), levels=np.linspace(Omin, Omax, levels), extend=extend, cmap=cmap)
    ax.gridlines(crs=ccrs.PlateCarree())
    if cbr:
        cb1 = plt.colorbar(h1, ax=ax, fraction=0.035, aspect=19, orientation="horizontal", pad=0.1)
        cb1.set_label(r"ODF$(\theta,\phi)$ - mean(ODF) (i.e. isotropic contr. removed)")
    if show:
        print(np.mean(ODF))
        plt.show()
    return h1


def fabric_to_hor_rot(f1, f2, f5):
    """Skip the overly complicated calculations, just go horizontal.

    Essentially assumes that one of the E_i's is vertical.
    Result in degrees.
    """
    A_xy = np.zeros((len(f1), 2, 2))
    A_xy[:, 0, 0] = f1
    A_xy[:, 0, 1] = f5
    A_xy[:, 1, 0] = f5
    A_xy[:, 1, 1] = 1.0 - f1 - f2
    A_xy[np.isnan(A_xy)] = 1.0
    A_xy[np.isinf(A_xy)] = 1.0
    E, rot_mat = np.linalg.eig(A_xy)

    norm = (A_xy[:, 0, 0] + A_xy[:, 1, 1]) / (E[:, 0] + E[:, 1])
    E[:, 0] = E[:, 0] * norm
    E[:, 1] = E[:, 1] * norm
    φ = np.arccos(-rot_mat[:, 0, 0]) * 180.0 / np.pi
    φ[rot_mat[:, 0, 1] < 0.0] = -φ[rot_mat[:, 0, 1] < 0.0]
    φ[φ < -90.0] = φ[φ < -90.0] + 180.0
    φ[φ > 90.0] = φ[φ > 90.0] - 180.0

    ## These are really unnecessary
    # φ[φ < -45.0] = φ[φ < -45.0] + 90.0
    # φ[φ > 45.0] = φ[φ > 45.0] - 90.0

    return φ


def fabric_to_hor_eigs(f1, f2, f5, elmer_ordering=True, normalize=False):
    """Skip the overly complicated calculations, just go horizontal.

    Essentially assumes that one of the E_i's is vertical.
    Result in degrees.

    Parameters
    ----------
    f1: np.array
        The xx fabric (f1 for Elmer's world).
    f2: np.array
        If elmer_ordering, f2 is the zz fabric (f2 in Elmer; yy calculated). Else, the yy fabric.
    f5: np.array
        The xy fabric (f5 in Elmer's world).
    elmer_ordering: bool, optional
        Is f2 zz (true) or yy (false). Default zz (so can take direct output from Elmer).
    normalize: bool, optional
        If false, the horizontal output is going to sum to 1. Otherwise, sums to 1 - f_zz
    """
    A_xy = np.zeros((len(f1), 2, 2))
    A_xy[:, 0, 0] = f1
    A_xy[:, 0, 1] = f5
    A_xy[:, 1, 0] = f5
    A_xy[:, 1, 1] = 1.0 - f1 - f2
    A_xy[np.isnan(A_xy)] = 1.0
    A_xy[np.isinf(A_xy)] = 1.0
    E, rot_mat = np.linalg.eig(A_xy)

    norm = (A_xy[:, 0, 0] + A_xy[:, 1, 1]) / (E[:, 0] + E[:, 1])
    E[:, 0] = E[:, 0] * norm
    E[:, 1] = E[:, 1] * norm

    φ = np.arccos(-rot_mat[:, 0, 0]) * 180.0 / np.pi
    φ[rot_mat[:, 0, 1] < 0.0] = -φ[rot_mat[:, 0, 1] < 0.0]
    φ[φ < -90.0] = φ[φ < -90.0] + 180.0
    φ[φ > 90.0] = φ[φ > 90.0] - 180.0

    ## These are really unnecessary
    # φ[φ < -45.0] = φ[φ < -45.0] + 90.0
    # φ[φ > 45.0] = φ[φ > 45.0] - 90.0

    if normalize:
        if elmer_ordering:
            A_zz = f2
        else:
            A_zz = 1.0 - f1 - f2
        target_totals = 1.0 - A_zz
        E[:, 0] = target_totals * E[:, 0]
        E[:, 1] = target_totals * E[:, 1]
    return φ, E


def fabric_to_ver_rot(f1, f2, f3, quadrant=True):
    """Skip the overly complicated calculations, just go horizontal.

    Essentially assumes that one of the E_i's is vertical.
    Result in degrees.
    """
    A_xy = np.zeros((len(f1), 2, 2))
    A_xy[:, 0, 0] = f1
    A_xy[:, 0, 1] = f3
    A_xy[:, 1, 0] = f3
    A_xy[:, 1, 1] = f2
    A_xy[np.isnan(A_xy)] = 1.0
    A_xy[np.isinf(A_xy)] = 1.0
    E, rot_mat = np.linalg.eig(A_xy)

    # Rotation ccw from the x axis
    φ = np.arctan2(rot_mat[:, 1, 0], rot_mat[:, 0, 0]) * 180.0 / np.pi

    # Do this mod 180 since axes are bidirectional
    φ[φ < -90.0] = φ[φ < -90.0] + 180.0
    φ[φ > 90.0] = φ[φ > 90.0] - 180.0

    if quadrant:
        # We are going to take the one closest to vertical (we can do this because of orthogonality)
        φ[φ < -45.0] = φ[φ < -45.0] + 90.0
        φ[φ > 45.0] = φ[φ > 45.0] - 90.0

    return φ


def quiver(ax, vtu, s=75, x=np.arange(-83333, 100000, 33333.0), y=np.arange(100, 2000, 400.0), xshift=lambda x: x / 1000.0 + 100.0, scale=1.0, width=None, zorder=2):
    X, Y = np.meshgrid(x, y)
    dat = vtu.get_pts_2d(["eigenv 1", "eigenv 2", "eigenv 3", "eigenv 4", "fabric 1", "fabric 2", "fabric 3", "fabric 5"], X.flatten(), Y.flatten())
    return quiver_dict(ax, dat, s=s, X=xshift(X), Y=Y, scale=scale, width=width, zorder=zorder)


def quiver_dict(ax, dat, s=75, X=None, Y=None, inx=True, scale=1.0, width=None, zorder=2):
    hang = fabric_to_ver_rot(dat["fabric 1"], dat["fabric 2"], dat["fabric 5"])
    if inx:
        ang = fabric_to_ver_rot(dat["fabric 1"], dat["fabric 2"], dat["fabric 3"])
        ang[ang < -45] = ang[ang < -45] + 90.0
        ang[ang > 45] = ang[ang > 45] - 90.0
        u = np.sin(ang / 180.0 * np.pi) * dat["eigenv 3"]
        v = np.cos(ang / 180.0 * np.pi) * dat["eigenv 3"]
        units = "width"
    else:
        u = np.zeros_like(dat["eigenv 3"])
        v = np.ones_like(dat["eigenv 3"])
        units = "height"
        ang = np.zeros_like(dat["fabric 1"])
    print(np.max(np.abs(hang)), np.max(np.abs(ang)))

    singlemax = dumb_woodcock(dat)
    vertical = dat["fabric 2"] > (1.0 - dat["fabric 2"] - dat["fabric 1"])
    vertical_sm = np.logical_and(vertical, singlemax)
    hor_sm = np.logical_and(~vertical, singlemax)
    quiv = ax.quiver(X.flatten()[vertical_sm], Y.flatten()[vertical_sm], u[vertical_sm], v[vertical_sm], units=units, scale=scale, width=width, zorder=zorder)
    if np.any(hor_sm):
        hq = [ax.plot(X.flatten()[hor_sm], Y.flatten()[hor_sm], marker=".", markersize=2, linestyle="none", color="0.4", zorder=zorder, label="Single max. partly into page")]
    else:
        hq = []

    planlabel = False
    ooplabel = False
    other_pts = []
    for i in range(np.sum(~singlemax)):
        t = MarkerStyle(marker=HOOP)
        t._transform = t.get_transform().rotate_deg(-ang[~singlemax][i])
        if np.isnan(dat["fabric 1"].flatten()[~singlemax][i]):
            continue
        if np.abs(hang.flatten()[~singlemax][i]) > 1:
            color = "0.4"
            if not ooplabel:
                label = "Vert. girdle, normal out of x-z"
                ooplabel = True
            else:
                label = None
        else:
            color = "k"
            if not planlabel:
                label = "Vert. girdle, normal in x-z"
                planlabel = True
            else:
                label = None

        other_pts.append(ax.scatter(X.flatten()[~singlemax][i], Y.flatten()[~singlemax][i], marker=t, s=s, linewidth=0.5, c="none", zorder=zorder, edgecolors=color, label=label))
    return [quiv] + hq + other_pts


def plot_axes(ax, geo=ccrs.Geodetic()):
    cax = "tab:red"
    ax.plot([0, -90, 0], [0, 0, 90], marker=".", ms=2, c=cax, transform=geo, linestyle="none")
    ax.text(0, 0, r"$x$", fontsize=8, c=cax, transform=geo, ha="left", va="center")  # x axis
    ax.text(-90, 0, r"$-y$", fontsize=8, c=cax, transform=geo, ha="right", va="center")  # y axis
    ax.text(0, 90, r"$z$", fontsize=8, c=cax, transform=geo, ha="left", va="top")  # z axis


if __name__ == "__main__":
    # Demos
    # a2 = [[0, 0, 0], [0, 0, 0], [0, 0, 1]]  # max in z
    # a2 = [[1, 0, 0], [0, 0, 0], [0, 0, 0]]  # max in x
    # a2 = [[0, 0, 0], [0, 1, 0], [0, 0, 0]]  # max in y
    a2 = [[1 / np.sqrt(2), 0, 0], [0, 1 / np.sqrt(2), 0], [0, 0, 0]]  # x-y girdle
    # a2 = [[1 / np.sqrt(2), 0, 0], [0, 0, 0], [0, 0, 1 / np.sqrt(2)]]  # x-z girdle
    # a2 = [[1. / 3., 0, 0], [0, 1. / 3., 0], [0, 0, 1. / 3.]]  # iso

    a2plot(np.matrix(a2))
