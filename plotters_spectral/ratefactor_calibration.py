#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <david.lilien@umanitoba.ca>, 2022
#
# Distributed under terms of the GNU GPL3.0 license.


import numpy as np
from scipy.interpolate import interp1d
import pickle
from progress.bar import Bar
import matplotlib.pyplot as plt

from specfabpy import specfabpy as sf
from ratefactor_calibration_dj import setup_model

CALIBRATE = False

FS = 10
BFS = 14
cm = plt.get_cmap("tab20")


def c2k(degc):
    return degc + 273.15  # deg. C to deg. K


class PureShear:
    def __init__(self, t_e, r=0, ax="z"):
        self.t_e = float(t_e)  # e-folding time scale for for parcel height reduction
        if ax == "z":
            self.Fpow = [(1.0 + r) / 2.0, (1.0 - r) / 2.0, -1]  # r=0 => unconfined

    def lam(self, t):
        return np.exp(t / self.t_e)  # lambda(t)

    def F(self, t):
        return np.diag(np.power(self.lam(t), self.Fpow))  # Deformation tensor

    def strain(self, t):
        return 0.5 * (self.F(t) + np.transpose(self.F(t))) - np.diag([1.0, 1, 1])  # Strain tensor

    def strainzz2time(self, strain_zz):
        return -self.t_e * np.log(strain_zz + 1)  # time it takes to reach "strain_zz" strain with timescale t_e.

    def D(self):
        return 1 / self.t_e * np.diag(self.Fpow)  # Strain rate. Note that F is constructed such that W and eps are time-independant.

    def W(self):
        return np.diag([0.0, 0, 0])  # Spin


def f_lam(D, T, m_lam=1.26e-3, b_lam=0.206):
    eps_E = np.sqrt(0.5 * np.einsum("ij,ji", D, D))  # sqrt(0.5 * D:D)
    return max(eps_E * (m_lam * T + b_lam), 0.0)


L_list = [8]
L = 8
lm, nlm_len = sf.init(L)
PROFILE = "EDC"

Nt = 200  # Number of integration steps

with open("DOME_C.dump", "rb") as fin:
    (fab, temp, aux) = pickle.load(fin)
n22_0 = -0.01081568
n20_0 = 0.10218687
strain_zz_stop = -0.96  # stop early, fabric only measure at EDC until roughly this point.

# Geometry
z_bed = -aux["H"]  # depth of bed
z0 = abs(fab["z"][0])
z1 = aux["H"]
H = z1 - z0  # ice column height as seen by parcel model

# Fabric
lami_meas, z_meas = fab["eig"], fab["z"]

# Temperature interpolant
T, zT = temp["T"], temp["z"]
f_T = interp1d(zT, T)

# Accumulation rate
# b = aux['b']
# Use the mean accumulation from EDC, not modern
b = 0.0153

t_e = H / b  # characteristic timescale
ps = PureShear(t_e, r=0.0, ax="z")  # r=0 => unconfined

# Constants
tend = ps.strainzz2time(strain_zz_stop)
timevec = np.linspace(0, tend, Nt)
dt = timevec[1] - timevec[0]  # given t_e, if we want Nt time steps, this is the step size needed
# Assumes depth-constant (unconfined) vertical compression: the Nye's classical dome model.
strain_zz = np.array([ps.strain(t)[-1, -1] for ii, t in enumerate(timevec)])  # for the constant time-step size, these are the vertical parcel strains as a function of time
z_sf = H * (strain_zz + 1) - H  # Vertical parcel strains correspond to these depths
z_sf += -z0  # effective offset of parcel model
T_sf = f_T(z_sf)


def run_model(A, Q, m_lam, b_lam, m_iota=0.0, c_iota=1.0, linear_mig=False):
    # Initialize arrays
    nlm = np.zeros((Nt, nlm_len), dtype=np.complex64)  # array of expansion coefficients
    nlm[0, 0] = 1.0 / np.sqrt(4 * np.pi)  # normalized distribution
    nlm[0, 3] = n20_0  # initial single max strength (fitted to match eigen value profile)
    nlm[0, 1] = n22_0  # initial single max strength (fitted to match eigen value profile)
    nlm[0, 5] = n22_0  # initial single max strength (fitted to match eigen value profile)
    lami_sf = np.zeros((Nt, 3))  # array of modelled eigen values
    lami_sf[0, :] = np.diag(sf.a2(nlm[0, :]))  # a2 already diagional (in eigenbasis) for this mode of deformation, so simply extract the diagonal values
    gam, lam = np.zeros((Nt)), np.zeros((Nt))

    # Euler integration of fabric model
    D, W = ps.D(), ps.W()  # strain-rate and spin tensors for mode of deformation
    if linear_mig:
        gam[0] = f_lam(D, T_sf[0], Q, A)
    else:
        gam[0] = sf.Gamma0(D, c2k(T_sf[0]), A, Q)

    lam[0] = f_lam(D, T_sf[0], m_lam=m_lam, b_lam=b_lam)

    with Bar("dt={:f}yr, Nt={:d} :: L={:d} (nlm_len={:d}) ::".format(dt, Nt, L, nlm_len), max=Nt - 1, fill="#", suffix="%(percent).1f%% - %(eta)ds") as bar:
        for nn in np.arange(1, Nt):
            nlm_prev = nlm[nn - 1, :]
            M_LATROT = (m_iota * T_sf[nn] + c_iota) * sf.M_LROT(nlm_prev, D, W, 1, 0)  # Lattice rotation operator (nlm_len x nlm_len matrix)
            S = D.copy()  # S remains coaxial with D for this mode of deformation, so we need not calculate S from the bulk flow law (which in turn requires calculating the enhancement factors from the modelled fabric)
            if linear_mig:
                gam[nn] = f_lam(D, T_sf[nn], Q, A)
            else:
                gam[nn] = sf.Gamma0(D, c2k(T_sf[nn]), A, Q)  # DDRX decary rate magnitude
            lam[nn] = f_lam(D, T_sf[nn], m_lam=m_lam, b_lam=b_lam)

            M_DDRX = gam[nn] * sf.M_DDRX(nlm_prev, S)  # DDRX operator (nlm_len x nlm_len matrix)
            M_CDRX = lam[nn] * sf.M_CDRX(nlm_prev)  # CDRX operator (nlm_len x nlm_len matrix)
            M_REG = sf.M_REG(nlm_prev, D)  # Regularization operator (nlm_len x nlm_len matrix)
            M = M_LATROT + M_DDRX + M_CDRX + M_REG  # Total fabric evolution operator (matrix)
            nlm[nn, :] = nlm_prev + dt * np.matmul(M, nlm_prev)  # Forward Euler step
            nlm[nn, :] = sf.apply_bounds(nlm[nn, :])  # Apply spectral bounds if needed
            lami_sf[nn, :] = np.diag(sf.a2(nlm[nn, :]))  # a2 already diagional (in eigen basis) for this mode of deformation, so simply extract the diagonal values
            bar.next()
    return z_sf, lami_sf, nlm, gam, lam


def misfit_function(z, lami, shallow=False):
    """Define the misfit just off the largest eigenvalue

    By symmetry, we always have l1==l2, so this is sufficient in this case.
    """
    lam3int = interp1d(z, lami[:, 2])
    if shallow:
        mask = z_meas > -3000.0
    else:
        mask = ~np.isnan(z_meas)
    return np.sqrt(np.mean(((lam3int(z_meas[mask]) - lami_meas[mask, 2])) ** 2.0)), np.sqrt(np.mean(((lam3int(z_meas[mask]) - lami_meas[mask, 2]) / lam3int(z_meas[mask])) ** 2.0)) * 100.0


def plot_single_result(ax, lami_sf, z_sf, name, ls="solid"):
    ax.plot(lami_sf[:, 0], -z_sf, c="C2", ls=ls, lw=2, clip_on=False)
    ax.plot(lami_sf[:, 1], -z_sf, c="C1", ls=ls, lw=2, clip_on=False)
    ax.plot(lami_sf[:, 2], -z_sf, c="C0", ls=ls, lw=2, clip_on=False)
    ax.plot([-1, -1], [1, 1], label=name, c="k", ls=ls, lw=2, clip_on=False)


# Q = 3.36e4
# A = 1.91e7
Q = 0.176
A = 6.09
m_lam = 1.26e-3
b_lam = 0.206
m_iota = 0.026
c_iota = 1.95

z_iota, lami_iota, nlm, gam, lam = run_model(A, Q, m_lam, b_lam, m_iota, c_iota, linear_mig=True)
print("Iota misfit is {:f} ({:f}%)".format(*misfit_function(z_iota, lami_iota)))

Q = 3.36e4
m_lam = 1.26e-3
b_lam = 0.206
A = 1.91e7
m_iota = 0.0
c_iota = 1.0

z_exp, lami_exp, nlm, gam, lam = run_model(A, Q, m_lam, b_lam, m_iota, c_iota)
print("Experimental param. misfit is {:f} ({:f}%)".format(*misfit_function(z_exp, lami_exp)))

Q = 3.36e4
m_lam = 0.0
b_lam = 0.0
A = 4.4e7
m_iota = 0.0
c_iota = 1.0

z_cal, lami_cal, nlm, gam, lam = run_model(A, Q, m_lam, b_lam, m_iota, c_iota)
print("Calibrated param. misfit is {:f} ({:f}%)".format(*misfit_function(z_cal, lami_cal)))
print("Calibrated param. misfit shallow is {:f} ({:f}%)".format(*misfit_function(z_cal, lami_cal, shallow=True)))

run_model_dj_02 = setup_model(0.2)
run_model_dj_04 = setup_model(0.4)

m_lam = 1.26e-3
if CALIBRATE:
    m_lam = 1.26e-3
    As = np.arange(0.0, 1.0e8, 1.0e6)
    bs = np.arange(0.0, 0.1, 0.001)
    misfits = np.zeros((As.shape[0], bs.shape[0]))
    for i, A in enumerate(As):
        for j, b_lam in enumerate(bs):
            z_shallow, lami_shallow, nlm, gam, lam = run_model(A, Q, m_lam, b_lam, m_iota, c_iota)
            misfits[i, j] = misfit_function(z_shallow, lami_shallow, shallow=True)[0]

    plt.figure()
    plt.title("shallow")
    cmm = plt.imshow(np.flipud(misfits), extent=(np.min(bs), np.max(bs), np.min(As), np.max(As)))
    plt.colorbar(cmm)
    plt.xlabel(r"$b_\Lambda$")
    plt.ylabel(r"$A_\Gamma$")
    plt.axis("auto")
    plt.savefig("../plots/raster_calibration_misfit_shallow.png")

    i, j = np.unravel_index(np.nanargmin(misfits, axis=None), misfits.shape)
    A_shallow = As[i]
    b_lam_shallow = bs[j]

    Q = 3.36e4
    m_lam = 1.26e-3
    b_lam = 4.0e-2
    A = 4.2e7
    m_iota = 0.0
    c_iota = 1.0

    As = np.arange(0.0, 1.0e8, 1.0e6)
    bs = np.arange(0.0, 0.1, 0.001)
    misfits = np.zeros((As.shape[0], bs.shape[0]))
    for i, A in enumerate(As):
        for j, b_lam in enumerate(bs):
            z_dj_02, lami_dj_02, nlm, gam, lam = run_model_dj_02(A, Q, m_lam, b_lam, m_iota, c_iota)
            misfits[i, j] = misfit_function(z_dj_02, lami_dj_02, shallow=True)[0]

    plt.figure()
    plt.title("D--J")
    cmm = plt.imshow(np.flipud(misfits), extent=(np.min(bs), np.max(bs), np.min(As), np.max(As)))
    plt.colorbar(cmm)
    plt.xlabel(r"$b_\Lambda$")
    plt.ylabel(r"$A_\Gamma$")
    plt.axis("auto")
    plt.savefig("../plots/raster_calibration_misfit_dj02.png")

    i, j = np.unravel_index(np.nanargmin(misfits, axis=None), misfits.shape)
    A_dj_02 = As[i]
    b_lam_dj_02 = bs[j]

    As = np.arange(0.0, 1.0e8, 1.0e6)
    bs = np.arange(0.0, 0.1, 0.001)
    misfits = np.zeros((As.shape[0], bs.shape[0]))
    for i, A in enumerate(As):
        for j, b_lam in enumerate(bs):
            z_dj_04, lami_dj_04, nlm, gam, lam = run_model_dj_04(A, Q, m_lam, b_lam, m_iota, c_iota)
            misfits[i, j] = misfit_function(z_dj_04, lami_dj_04, shallow=True)[0]

    plt.figure()
    plt.title("D--J")
    cmm = plt.imshow(np.flipud(misfits), extent=(np.min(bs), np.max(bs), np.min(As), np.max(As)))
    plt.colorbar(cmm)
    plt.xlabel(r"$b_\Lambda$")
    plt.ylabel(r"$A_\Gamma$")
    plt.axis("auto")
    plt.savefig("../plots/raster_calibration_misfit_dj04.png")

    i, j = np.unravel_index(np.nanargmin(misfits, axis=None), misfits.shape)
    A_dj_04 = As[i]
    b_lam_dj_04 = bs[j]
else:
    Q = 3.36e4
    m_lam = 1.26e-3
    b_lam = 4.0e-2
    A = 4.2e7
    m_iota = 0.0
    c_iota = 1.0

    A_dj_02, b_lam_dj_02 = 4.0e7, 3.2e-3
    z_dj_02, lami_dj_02, _, _, _ = run_model_dj_04(A_dj_02, Q, m_lam, b_lam_dj_02, m_iota, c_iota)
    A_dj_04, b_lam_dj_04 = 4.7e7, 5.2e-3
    z_dj_04, lami_dj_04, _, _, _ = run_model_dj_04(A_dj_04, Q, m_lam, b_lam_dj_04, m_iota, c_iota)
    A_shallow, b_lam_shallow = 4.6e7, 1.8e-3
    z_shallow, lami_shallow, nlm, gam, lam = run_model(A_shallow, Q, m_lam, b_lam_shallow, m_iota, c_iota)

    print("For the shallow version, A={:e}, b={:e}".format(A_shallow, b_lam_shallow))
    print("Shallow param. misfit is {:f} ({:f}%)".format(*misfit_function(z_shallow, lami_shallow, shallow=True)))
    print("For the DJ kink=0.2 version, A={:e}, b={:e}".format(A_dj_02, b_lam_dj_02))
    print("DJ 0.2 param. misfit is {:f} ({:f}%)".format(*misfit_function(z_dj_02, lami_dj_02, shallow=True)))
    print("For the DJ kink=0.4 version, A={:e}, b={:e}".format(A_dj_04, b_lam_dj_04))
    print("DJ 0.4 param. misfit is {:f} ({:f}%)".format(*misfit_function(z_dj_04, lami_dj_04, shallow=True)))

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 4.0), sharey=True, gridspec_kw=dict(left=0.12, right=0.98, top=0.98, bottom=0.13, width_ratios=[1, 0.5]))
ms = 5
ax1.plot(lami_meas[:, 0], -z_meas, marker="s", markeredgecolor="none", markerfacecolor=cm(5), linestyle="none", markersize=ms)
ax1.plot(lami_meas[:, 1], -z_meas, marker="s", markeredgecolor="none", markerfacecolor=cm(3), linestyle="none", markersize=ms)
ax1.plot(lami_meas[:, 2], -z_meas, marker="s", markeredgecolor="none", markerfacecolor=cm(1), linestyle="none", markersize=ms)
ax1.plot(-1, -1, marker="s", markeredgecolor="none", markerfacecolor="k", linestyle="none", markersize=ms, label="EDC data")

ax1.set_xlabel(r"$\lambda_i$", fontsize=FS)
ax1.set_xlim([0, 1])

ax1.set_ylabel("Depth (m)", fontsize=FS)
ax1.set_ylim([-z_bed, 0])
ax1.set_yticks(np.flipud(np.arange(0, H, 500)))
ax1.set_yticks(np.flipud(np.arange(0, H, 250)), minor=True)

ax1.set_xticks([0.0, 1.0 / 3.0, 2.0 / 3.0, 1])
ax1.set_xticklabels(["0", r"$\frac{1}{3}$", r"$\frac{2}{3}$", "1"])

plot_single_result(ax1, lami_exp, z_exp, "Lab", ls="solid")
plot_single_result(ax1, lami_cal, z_cal, "Ice core", ls="dashed")
plot_single_result(ax1, lami_iota, z_iota, "R2021", ls="dotted")

ax1.tick_params(axis="both", which="major", labelsize=FS)
ax1.legend(loc="best", fontsize=FS, frameon=False)

ax2.set_xlabel(r"Temperature ($^\circ$C)")
ax2.plot(T, -zT, color="k")
fig.savefig("../plots/figure_2.pdf", dpi=300)

fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(5.5, 4.0), sharey=True, gridspec_kw=dict(left=0.12, right=0.98, top=0.98, bottom=0.13, width_ratios=[1, 0.5]))
ms = 5
ax1.plot(lami_meas[:, 0], -z_meas, marker="s", markeredgecolor="none", markerfacecolor=cm(5), linestyle="none", markersize=ms)
ax1.plot(lami_meas[:, 1], -z_meas, marker="s", markeredgecolor="none", markerfacecolor=cm(3), linestyle="none", markersize=ms)
ax1.plot(lami_meas[:, 2], -z_meas, marker="s", markeredgecolor="none", markerfacecolor=cm(1), linestyle="none", markersize=ms)
ax1.plot(-1, -1, marker="s", markeredgecolor="none", markerfacecolor="k", linestyle="none", markersize=ms, label="EDC data")

ax1.set_xlabel(r"$\lambda_i$", fontsize=FS)
ax1.set_xlim([0, 1])

ax1.set_ylabel("Depth (m)", fontsize=FS)
ax1.set_ylim([-z_bed, 0])
ax1.set_yticks(np.flipud(np.arange(0, H, 500)))
ax1.set_yticks(np.flipud(np.arange(0, H, 250)), minor=True)

ax1.set_xticks([0.0, 1.0 / 3.0, 2.0 / 3.0, 1])
ax1.set_xticklabels(["0", r"$\frac{1}{3}$", r"$\frac{2}{3}$", "1"])

# plot_single_result(ax1, lami_exp, z_exp, 'Lab', ls='solid')
plot_single_result(ax1, lami_cal, z_cal, "Full\ndepth", ls="dashed")
# plot_single_result(ax1, lami_iota, z_iota, 'R2021', ls='dotted')
plot_single_result(ax1, lami_shallow, z_shallow, "Limited\ndepth", ls="solid")
plot_single_result(ax1, lami_dj_02, z_dj_02, "Kink=0.2", ls="dotted")
plot_single_result(ax1, lami_dj_04, z_dj_04, "Kink=0.4", ls="dashdot")

ax1.tick_params(axis="both", which="major", labelsize=FS)
ax1.legend(loc="best", fontsize=FS, frameon=False)

ax2.set_xlabel(r"Temperature ($^\circ$C)")
ax2.plot(T, -zT, color="k")
fig.savefig("../plots/supplementary_figure_1.pdf", dpi=300)
