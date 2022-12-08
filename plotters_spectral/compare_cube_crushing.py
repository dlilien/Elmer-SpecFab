#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Plot up elmer cube experiments
"""
import numpy as np
import pyvista
import matplotlib.pyplot as plt
import matplotlib.ticker as ticker
import matplotlib.gridspec as gridspec
import fabricplotlib
from specfabpy import specfabpy as sf
import netCDF4
import cartopy.crs as ccrs

from load_experimental import load_pickle

experiments = load_pickle(False)


plt.rcParams["ytick.labelsize"] = 8
plt.rcParams["xtick.labelsize"] = 8
plt.rcParams["axes.labelsize"] = 10
plt.rcParams["axes.titlesize"] = 12

rot = -45
inclination = 70
prj = ccrs.Orthographic(rot, 90 - inclination)

LCap_elmer = 6
landm, nlm_len_elmer = sf.init(LCap_elmer)
nlm_len_elmer = nlm_len_elmer * 2
lvals = landm[0, :]
mvals = landm[1, :]

markerdict = {"Qi": "d", "Hunter": "o", "Fan": "X"}

cms = [None, None, plt.get_cmap("Reds"), None, plt.get_cmap("Blues"), None, plt.get_cmap("Greens"), None, plt.get_cmap("Oranges")]

fabplot_args = {"show": False, "cbr": False, "levels": 21, "Omin": 0.0, "Omax": 0.4, "latres": 100, "lonres": 100, "cmap": "Greys"}

FNAMES = ["DDRX", "FULL", "LATROT"]
rc_types = ["stress", "strain"]
UGRADS = {"CC_ZX": "Confined compression", "SS_XZ": "Simple shear", "UC_ZZ": "Unconfined compression", "UE_ZZ": "Uniform extension"}
temps = ["-30", "-5", "0"]

DIMS = {"CC_ZX": 2, "SS_XZ": 2, "UC_ZZ": 3, "UE_ZZ": 3}
rc_linestyles = ["dashed", "dotted"]
rc_linestyle_dicts = [{"linestyle": "none", "marker": "x", "markersize": 5}, {"linestyle": "dotted"}]

N_STEPS = 501
TIMES = np.arange(0.0, N_STEPS)


def plot_eigs(ax, x, a2, marker="o", markersize=5, linestyle=None, **kwargs):
    eigs = np.linalg.eigvals(a2)
    eigs = np.sort(eigs)

    # We can have problems with negative values from finite precision
    if eigs[0] < 0.0:
        eigs[0] = 0.0
    if eigs[1] < 0.0:
        eigs[1] = 0.0
    if eigs[2] > 1.0:
        eigs[2] = 1.0
    eigs = eigs / np.sum(np.sqrt(eigs**2.0))
    if linestyle is not None:
        ax.axhline(eigs[0], color="C2", linestyle=linestyle, **kwargs)
        ax.axhline(eigs[1], color="C1", linestyle=linestyle, **kwargs)
        ax.axhline(eigs[2], color="C0", linestyle=linestyle, **kwargs)
    else:
        ax.plot(x, eigs[0], color="C2", marker=marker, markersize=markersize, linestyle=None, **kwargs)
        ax.plot(x, eigs[1], color="C1", marker=marker, markersize=markersize, linestyle=None, **kwargs)
        ax.plot(x, eigs[2], color="C0", marker=marker, markersize=markersize, linestyle=None, **kwargs)


def unified_plot(full_fabs):
    fig = plt.figure(figsize=(7, 4.25))
    gs = gridspec.GridSpec(3, 12, wspace=0.0, hspace=0.1, left=0.055, right=0.985, bottom=0.15, top=0.96, width_ratios=[1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 1, 0.1, 0.06, 0.35], height_ratios=[1, 1, 0.75])
    top_ball_axes = [fig.add_subplot(gs[0, i * 2], projection=prj) for i in range(5)]
    bot_ball_axes = [fig.add_subplot(gs[1, i * 2], projection=prj) for i in range(5)]

    gs0 = gridspec.GridSpecFromSubplotSpec(1, 4, subplot_spec=gs[2, :])

    marker_axes = [fig.add_subplot(gs0[0, 0])]
    for i in range(1, 4):
        marker_axes.append(fig.add_subplot(gs0[0, i]))  # , sharex=marker_axes[0]))

    cax = fig.add_subplot(gs[:-1, 10])

    # marker_ax = fig.add_subplot(gs[:-1, -1])
    temp = "-5"
    fabricplotlib.nlm_plot(LCap_elmer, full_fabs[temp]["stress"]["LATROT"]["SS_XZ"][-1, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp]["stress"]["LATROT"]["SS_XZ"][-1, nlm_len_elmer // 2 :], elmer=False, ax=top_ball_axes[0], **fabplot_args)
    cm = fabricplotlib.nlm_plot(
        LCap_elmer, full_fabs[temp]["stress"]["LATROT"]["UE_ZZ"][-1, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp]["stress"]["LATROT"]["UE_ZZ"][-1, nlm_len_elmer // 2 :], elmer=False, ax=bot_ball_axes[0], **fabplot_args
    )
    top_ball_axes[0].text(0.5, 1.1, "Lat. Rot.", ha="center", va="center", fontsize=10, transform=top_ball_axes[0].transAxes)

    for i, temp in enumerate(temps):
        fabricplotlib.nlm_plot(
            LCap_elmer, full_fabs[temp]["stress"]["FULL"]["SS_XZ"][-1, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp]["stress"]["FULL"]["SS_XZ"][-1, nlm_len_elmer // 2 :], elmer=False, ax=top_ball_axes[i + 1], **fabplot_args
        )
        cm = fabricplotlib.nlm_plot(
            LCap_elmer, full_fabs[temp]["stress"]["FULL"]["UE_ZZ"][-1, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp]["stress"]["FULL"]["UE_ZZ"][-1, nlm_len_elmer // 2 :], elmer=False, ax=bot_ball_axes[i + 1], **fabplot_args
        )
        top_ball_axes[i + 1].text(0.5, 1.1, temp + r"$^\circ$C", ha="center", va="center", fontsize=10, transform=top_ball_axes[i + 1].transAxes)

    temp = "-5"
    fabricplotlib.nlm_plot(LCap_elmer, full_fabs[temp]["stress"]["DDRX"]["SS_XZ"][-1, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp]["stress"]["DDRX"]["SS_XZ"][-1, nlm_len_elmer // 2 :], elmer=False, ax=top_ball_axes[4], **fabplot_args)
    cm = fabricplotlib.nlm_plot(
        LCap_elmer, full_fabs[temp]["stress"]["DDRX"]["UE_ZZ"][-1, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp]["stress"]["DDRX"]["UE_ZZ"][-1, nlm_len_elmer // 2 :], elmer=False, ax=bot_ball_axes[4], **fabplot_args
    )
    top_ball_axes[4].text(0.5, 1.1, "DDRX", ha="center", va="center", fontsize=10, transform=top_ball_axes[4].transAxes)

    for ax, letter in zip(top_ball_axes + bot_ball_axes, "abcdefghijkl"):
        fabricplotlib.plot_axes(ax)
        ax.text(0.05, 0.95, "{:s}.".format(letter), ha="center", va="center", fontsize=12, transform=ax.transAxes, fontweight="bold")

    top_ball_axes[0].text(-0.1, 0.5, "Simple shear", ha="center", va="center", fontsize=10, rotation=90, transform=top_ball_axes[0].transAxes)
    bot_ball_axes[0].text(-0.15, 0.5, "Uniform\nextension", ha="center", va="center", fontsize=10, rotation=90, transform=bot_ball_axes[0].transAxes)

    plt.colorbar(cm, cax=cax, orientation="vertical", label=r"ODF($\theta,\phi$)", ticks=[0, 0.1, 0.2, 0.3, 0.4, 0.5])

    for ugrad, marker_ax, letter in zip(UGRADS, marker_axes, "klmnop"):
        As = np.zeros((3, len(temps)))
        for i, temp in enumerate(temps):
            a2 = full_a2s[temp]["stress"]["FULL"][ugrad][-1, :, :]
            eigs = np.linalg.eigvals(a2)
            As[:, i] = np.sort(eigs)

        for i in range(3):
            marker_ax.plot([float(t) for t in temps], As[i, :], color="C{:d}".format(2 - i), marker="s", markersize=5, markeredgecolor="k")

        for expr in experiments:
            if expr.strain_type.upper() == ugrad[:2]:
                eigs = np.sort(np.linalg.eigvals(expr.a2))
                for i in range(3):
                    if (expr.strain_type == "ss" and expr.tot_strain > 2.0) or (expr.tot_strain > 0.35):
                        fc = "k"
                    else:
                        fc = "w"
                    marker_ax.plot(expr.temp, eigs[i], markeredgecolor="C{:d}".format(2 - i), marker=markerdict[expr.src], markersize=5, markerfacecolor=fc)
        marker_ax.set_xlim(-31, 0.5)
        marker_ax.set_ylim(0, 1)
        marker_ax.set_yticks([0, 1.0 / 3.0, 2.0 / 3.0, 1])
        marker_ax.set_yticklabels(["0", r"$\frac{1}{3}$", r"$\frac{2}{3}$", "1"])
        marker_ax.text(0.01, 0.98, r"$\bf{{{}.}}$ {:s}".format(letter, UGRADS[ugrad].replace(" ", "\n")), transform=marker_ax.transAxes, ha="left", va="top", fontsize=9)

    marker_axes[2].text(-0.1, -0.25, r"Temperature ($^\circ$C)", transform=marker_axes[2].transAxes, ha="center", va="top", fontsize=8)
    marker_axes[0].set_ylabel("Eigenvalues")

    marker_axes[0].plot(100, 100, marker="s", markersize=5, color="k", label="Model")
    marker_axes[0].plot(100, 100, marker=markerdict["Qi"], markersize=5, color="k", linestyle="none", label="Qi et al., 2019")
    marker_axes[0].plot(100, 100, marker=markerdict["Fan"], markersize=5, color="k", linestyle="none", label="Fan et al., 2020")
    marker_axes[0].plot(100, 100, marker=markerdict["Hunter"], markersize=5, color="k", linestyle="none", label="Hunter et al., 2022")

    marker_axes[0].legend(ncol=5, loc="upper left", fontsize=8, bbox_to_anchor=(-0.1, -0.45), frameon=False)

    fig.savefig("../plots/figure_9.png", dpi=300)


def eigs_notensor_plot(temp, full_a2s, spectral_a2s):
    fig = plt.figure(figsize=(7.047, 3.75))
    gs = gridspec.GridSpec(1, 4, wspace=0.2, hspace=0.0, left=0.055, right=0.99, bottom=0.11, top=0.88)
    axes = [fig.add_subplot(gs[0, i]) for i in range(4)]

    for ugrad, ax, letter in zip(UGRADS, axes, "abcdefg"):
        ax.set_title(r"$\mathbf{{{}.}}$".format(letter) + UGRADS[ugrad].replace(" ", "\n"), fontsize=12, loc="left")

        plot_eigs(ax, 1.3, full_a2ts[ugrad][-1, :, :], linestyle="dashed", zorder=0.1)

        # just do these manually since there is a lot of logic
        plot_eigs(ax, -0.3, spectral_a2s["FULL"]["L20"][ugrad][-1, :, :], marker="o", markeredgecolor="k")
        plot_eigs(ax, 0.0, full_a2s[temp]["stress"]["FULL"][ugrad][-1, :, :], marker="s", markeredgecolor="k")
        plot_eigs(ax, 0.3, full_a2s[temp]["strain"]["FULL"][ugrad][-1, :, :], marker="<", markeredgecolor="k")

        plot_eigs(ax, 0.7, spectral_a2s["LATROT"]["L20"][ugrad][-1, :, :], marker="o", markeredgecolor="k")
        plot_eigs(ax, 1.0, full_a2s[temp]["stress"]["LATROT"][ugrad][-1, :, :], marker="s", markeredgecolor="k")
        # Since DDRX is off, stress=strain DDRX and we only ran one
        plot_eigs(ax, 1.3, full_a2s[temp]["stress"]["LATROT"][ugrad][-1, :, :], marker="<", markeredgecolor="k")

        plot_eigs(ax, 1.7, spectral_a2s["DDRX"]["L20"][ugrad][-1, :, :], marker="o", markeredgecolor="k")
        plot_eigs(ax, 2.0, full_a2s[temp]["stress"]["DDRX"][ugrad][-1, :, :], marker="s", markeredgecolor="k")
        plot_eigs(ax, 2.3, full_a2s[temp]["strain"]["DDRX"][ugrad][-1, :, :], marker="<", markeredgecolor="k")

        for exp in experiments:
            if exp.temp > -8 and exp.temp < -3 and exp.strain_type.upper() == ugrad[:2]:
                plot_eigs(ax, 1.3, exp.a2, linestyle="solid", zorder=0.1)

                eigs = np.sort(np.linalg.eigvals(exp.a2))
                if exp.name == "D5-1":
                    ax.text(2.0, eigs[2] - 0.005, exp.name, fontsize=9, ha="center", va="top")
                else:
                    ax.text(2.0, eigs[2], exp.name, fontsize=9, ha="center", va="bottom")
                # for x in np.linspace(-0.5, 2.5, 20):
                #     plot_eigs(ax, 1.3, exp.a2, marker='x')

        ax.set_xlim(-0.5, 2.5)
        ax.set_xticks([0.5, 1.5])
        ax.set_xticks([0, 1, 2], minor=True)
        ax.grid(axis="x")
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
        ax.xaxis.set_tick_params(length=0, which="minor")
        ax.set_xticklabels(["Full", "Lat. Rot.", "DDRX"], minor=True)
        ax.set_ylim(0, 1)
        ax.set_yticks([0, 1.0 / 3.0, 2.0 / 3.0, 1])
        ax.set_yticklabels(["0", r"$\frac{1}{3}$", r"$\frac{2}{3}$", "1"])

    axes[0].plot(100, 100, marker="o", markersize=5, linestyle="none", color="k", label="SpecFab ($L=20$)")
    axes[0].plot(100, 100, marker="s", markersize=5, linestyle="none", color="k", label=r"Spectral, $\tau$ DDRX")
    axes[0].plot(100, 100, marker="<", markersize=5, linestyle="none", color="k", label=r"Spectral, $\dot{\varepsilon}$ DDRX")
    axes[0].plot([99, 100], [99, 100], linestyle="--", color="k", label=r"Tensorial (lat. rot.)")
    axes[0].plot([99, 100], [99, 100], linestyle="solid", color="k", label=r"Experimental")

    axes[0].legend(ncol=5, loc="upper left", fontsize=8, bbox_to_anchor=(-0.37, -0.05), frameon=False)

    axes[0].set_ylabel("Eigenvalues", fontsize=10)
    fig.savefig("../plots/figure_8.pdf", dpi=300)


def load_elmer_data(LCap_elmer):
    template = "../cube_crushing_{:s}rc/{:s}/{:s}_{:s}_{:s}_t{:04d}.vtu"

    full_fabs = {temp: {rc: {fname: {ugrad: np.zeros((N_STEPS, nlm_len_elmer)) for ugrad in UGRADS} for fname in FNAMES} for rc in rc_types} for temp in temps}
    full_a2s = {temp: {rc: {fname: {ugrad: np.zeros((N_STEPS, 3, 3)) for ugrad in UGRADS} for fname in FNAMES} for rc in rc_types} for temp in temps}
    full_a4s = {temp: {rc: {fname: {ugrad: np.zeros((N_STEPS, 3, 3, 3, 3)) for ugrad in UGRADS} for fname in FNAMES} for rc in rc_types} for temp in temps}

    # These cannot have processes, temp, or rc
    full_a2ts = {ugrad: np.zeros((N_STEPS, 3, 3)) for ugrad in UGRADS}

    for temp in temps:
        for rc in rc_types:
            for ugrad in UGRADS:
                for fname in FNAMES:
                    # just get the center
                    if DIMS[ugrad] == 2:
                        folder = "square"
                        try:
                            ss = pyvista.read(template.format(rc, folder, ugrad, fname, temp, 501))
                        except FileNotFoundError:
                            continue
                        ind_center = np.argmin((ss.points[:, 0] - 0.5) ** 2.0 + (ss.points[:, 1] - 0.5) ** 2.0)
                    else:
                        folder = "cube"
                        try:
                            ss = pyvista.read(template.format(rc, folder, ugrad, fname, temp, 501))
                        except FileNotFoundError:
                            continue
                        ind_center = np.argmin((ss.points[:, 0] - 0.5) ** 2.0 + (ss.points[:, 1] - 0.5) ** 2.0 + (ss.points[:, 2] - 0.5) ** 2.0)
                    for i in range(0, N_STEPS):
                        try:
                            vtu = pyvista.read(template.format(rc, folder, ugrad, fname, temp, i + 1))

                            # and the static plots
                            for j in range(nlm_len_elmer):
                                full_fabs[temp][rc][fname][ugrad][i, j] = vtu["{:s} {:d}".format("fabricr", j + 1)][ind_center]

                            full_a4s[temp][rc][fname][ugrad][i, :, :, :, :] = sf.a4(full_fabs[temp][rc][fname][ugrad][i, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp][rc][fname][ugrad][i, nlm_len_elmer // 2 :])

                            full_a2s[temp][rc][fname][ugrad][i, :, :] = sf.a2(full_fabs[temp][rc][fname][ugrad][i, : nlm_len_elmer // 2] + 1.0j * full_fabs[temp][rc][fname][ugrad][i, nlm_len_elmer // 2 :])

                            full_a2ts[ugrad][i, 0, 0] = vtu["{:s} {:d}".format("fabric", 1)][ind_center]
                            full_a2ts[ugrad][i, 1, 1] = vtu["{:s} {:d}".format("fabric", 2)][ind_center]
                            full_a2ts[ugrad][i, 2, 2] = 1.0 - full_a2ts[ugrad][i, 0, 0] - full_a2ts[ugrad][i, 1, 1]
                            full_a2ts[ugrad][i, 0, 1] = vtu["{:s} {:d}".format("fabric", 3)][ind_center]
                            full_a2ts[ugrad][i, 1, 0] = vtu["{:s} {:d}".format("fabric", 3)][ind_center]
                            full_a2ts[ugrad][i, 1, 2] = vtu["{:s} {:d}".format("fabric", 4)][ind_center]
                            full_a2ts[ugrad][i, 2, 1] = vtu["{:s} {:d}".format("fabric", 4)][ind_center]
                            full_a2ts[ugrad][i, 0, 2] = vtu["{:s} {:d}".format("fabric", 5)][ind_center]
                            full_a2ts[ugrad][i, 2, 0] = vtu["{:s} {:d}".format("fabric", 5)][ind_center]
                        except (FileNotFoundError, KeyError):
                            pass
    return full_a2ts, full_a2s, full_a4s, full_fabs


def load_spectral_data():
    # No RC or temp on these
    spectral_fab = {fname: {"L6": {}, "L8": {}, "L20": {}} for fname in FNAMES}
    spectral_a2s = {fname: {"L6": {}, "L8": {}, "L20": {}} for fname in FNAMES}
    spectral_a4s = {fname: {"L6": {}, "L8": {}, "L20": {}} for fname in FNAMES}
    for ugrad in UGRADS:
        for fname in FNAMES:
            for lc in spectral_fab[fname]:
                if DIMS[ugrad] == 2:
                    spectral_ds = netCDF4.Dataset("../cube_crushing_stressrc/specfab_0d/solutions_{:s}/{:s}_{:s}.nc".format(lc, fname, ugrad.lower().replace("z", "y")), mode="r")
                else:
                    spectral_ds = netCDF4.Dataset("../cube_crushing_stressrc/specfab_0d/solutions_{:s}/{:s}_{:s}.nc".format(lc, fname, ugrad.lower()), mode="r")

                def loadvar(field):
                    return np.array(spectral_ds.variables[field][:])

                # Model config
                spec_Nt = spectral_ds.getncattr("tsteps")

                # Fabric state
                spectral_fab[fname][lc][ugrad] = loadvar("c_re") + 1.0j * loadvar("c_im")

                # spectral_times = np.arange(0, spec_Nt * spec_dt, spec_dt)
                spectral_times = np.arange(0, spec_Nt)

                spectral_a2s[fname][lc][ugrad] = np.zeros((spectral_times.shape[0], 3, 3))
                spectral_a4s[fname][lc][ugrad] = np.zeros((spectral_times.shape[0], 3, 3, 3, 3))
                for i in range(spectral_a2s[fname][lc][ugrad].shape[0]):
                    spectral_a2s[fname][lc][ugrad][i, :, :] = sf.a2(spectral_fab[fname][lc][ugrad][i, :])
                    spectral_a4s[fname][lc][ugrad][i, :, :, :, :] = sf.a4(spectral_fab[fname][lc][ugrad][i, :])
    return spectral_times, spectral_a2s, spectral_a4s, spectral_fab


if __name__ == "__main__":
    spectral_times, spectral_a2s, spectral_a4s, spectral_fab = load_spectral_data()
    full_a2ts, full_a2s, full_a4s, full_fabs = load_elmer_data(LCap_elmer)
    spec_LCap = 6

    eigs_notensor_plot("-5", full_a2s, spectral_a2s)
    unified_plot(full_fabs)
