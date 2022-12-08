#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make figure 5 and supp. fig. 2
"""
import yaml
from glob import glob

import numpy as np
import matplotlib.pyplot as plt
from matplotlib.tri import Triangulation

from labellines import labelLines

import anisolib
import fabricplotlib
from modeltools.lib import fastvtulib


plt.rcParams["contour.negative_linestyle"] = "solid"
plt.rcParams.update({"mathtext.default": "regular"})


fs = 10
bfs = 14
EDC_DIST = 0.0

pres_fnum = 7


def domec_fabric():
    # We want to find the filename for the highest index where we have all 3
    fn_glob = "../domec_transect/domec_L6/domec_tube_full_aniso_calib_t????.vtu"
    files = glob(fn_glob)
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    fn = files[inds[-1]]
    fig, (ax1, ax2, ax3) = fabric(fn)
    ax3.set_xlabel("Distance from divide (km)", fontsize=fs)
    ax1.text(0.0, 1.01, "A", ha="center", va="bottom", fontsize=bfs, transform=ax1.transAxes)
    ax1.text(1.0, 1.01, "A'", ha="center", va="bottom", fontsize=bfs, transform=ax1.transAxes)

    vtu = fastvtulib.get_structured_vtu(fn)
    y_pts = 3000
    fabname = "tensorfabric"

    with open("pres_sites.yml") as fin:
        site_data = yaml.safe_load(fin)
    for site, letter in zip(site_data["sparse_sites"], "abcdefghijkl"):
        y = np.linspace(np.min(vtu.coords[:, 1]), np.max(vtu.coords[:, 1]), y_pts)
        x = np.ones((y_pts,)) * site_data["sites"][site]["dist"]
        test_p = vtu.get_pts_2d(["{:s} 1".format(fabname)], x, y)
        y_valid = y[~np.isnan(test_p["{:s} 1".format(fabname)])]
        y1 = np.max(y_valid)
        ax2.plot([site_data["sites"][site]["dist"] / 1000.0, site_data["sites"][site]["dist"] / 1000.0], [y1 - 2000.0, y1 - 1100.0], linestyle="dashed", color="firebrick", label=str(pres_fnum) + letter)
        ax2.plot([site_data["sites"][site]["dist"] / 1000.0, site_data["sites"][site]["dist"] / 1000.0], [y1 - 700.0, y1], linestyle="dashed", color="firebrick", label=str(pres_fnum) + letter)
        ax2.text(site_data["sites"][site]["dist"] / 1000.0, y1 - 1000.0, str(pres_fnum) + letter, fontsize=fs, color="firebrick", ha="center")
    labelLines(ax2.get_lines(), zorder=900, fontsize=fs, yoffsets=2500)
    fig.savefig("../plots/figure_5.png", dpi=300)


def divide_2p(background=""):
    xlim = (-15, 15)
    ylim = (0, 2025)
    overall_template = "../ideal_divide_spectral/divide_v2/{:s}_t{:s}.vtu"
    patterns = ["divide_ISO", "divide_spectral_guess", "divide_spectral_guess_GOLF"]
    rheos = ["Glen's flow law", "Full nonlinear\northotropic", "Nonlinear GOLF"]
    ind = 1e6
    for pattern in patterns:
        fn_glob = overall_template.format(pattern, "????")
        files = glob(fn_glob)
        inds = np.sort(np.array([int(fn[-8:-4]) for fn in files]))
        if inds[-1] < ind:
            ind = inds[-1]
    fig, (ax1, ax2) = plt.subplots(2, 1, figsize=(7.047, 5.5), sharex=True, sharey=True, gridspec_kw={"left": 0.1, "right": 0.985, "bottom": 0.16, "top": 0.99, "hspace": 0.05, "wspace": 0.0})

    for pattern, rheo, lc, ls in zip(patterns, rheos, ["C0", "C3", "C5"], ["solid", "solid", "solid"]):
        fn = overall_template.format(pattern, "{:04d}".format(ind))
        vtu = fastvtulib.get_structured_vtu(fn)

        tri = Triangulation(np.array([rc[0] / 1000.0 for rc in vtu.raw_coords[:, 0]]), np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt)

        levels = [5000, 10000, 20000, 30000, 40000, 50000, 75000]
        ctrs = ax1.tricontour(tri, vtu.rawdata_dict["age"], colors=lc, linestyles=ls, levels=levels, linewidths=1.0)
        fmt = {level: "{:d} ka".format(level // 1000) for level in levels}
        if rheo == rheos[0]:
            cls = plt.clabel(ctrs, fontsize=fs, inline=True, fmt=fmt, colors="k", manual=[(13.5, 1500), (13.5, 1250), (13.5, 800), (13.5, 500), (13.5, 340), (13.5, 240), (13.5, 100)])
            [txt.set_bbox(dict(facecolor="white", edgecolor="none", pad=0)) for txt in cls]

        # wlevels = [-0.01, -0.005, -0.001]
        # ax1.tricontour(tri, vtu.rawdata_dict['aiflow'][:, 1], colors=lc, linestyles='dashed', levels=wlevels, linewidths=0.5)

        ax1.plot([], [], label=rheo, color=lc, linestyle=ls, linewidth=1.0)

        levels = [-5e-5, -5e-6, 5e-6, 5e-5]
        ctrs = ax2.tricontour(tri, vtu.rawdata_dict["strainrate 4"], colors=lc, linestyles=ls, levels=levels, linewidths=1.0)
        plt.rcParams["font.size"] = "8"
        fmt = {level: r"$5\times10^{%d}$ a$^{-1}$" % np.log10(abs(level / 5)) for level in levels}
        for level in levels:
            if level < 0:
                fmt[level] = None
        if rheo == rheos[1]:
            cls = plt.clabel(ctrs, levels=[level for level in levels if level > 0], fontsize=fs, inline=True, fmt=fmt, colors="k", manual=[(12.5, 1800), (12.5, 900)])
            [txt.set_bbox(dict(facecolor="white", edgecolor="none", pad=2)) for txt in cls]
    ax1.legend(frameon=False, ncol=3, loc="upper left", fontsize=fs, bbox_to_anchor=(0.0, -1.25))

    for letter, ax in zip("abcdefghijkl", (ax1, ax2)):
        ax.text(0.01, 0.99, letter, ha="left", va="top", transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis="both", which="major", labelsize=fs)
        ax.set_ylim(*ylim)
        ax.set_xlim(*xlim)
        ax.set_ylabel("Elevation (m)", fontsize=fs)

    # for ax in (ax1):
    #     ax.tick_params(axis='x', which='both', labelbottom=False, labelleft=False)
    ax2.set_xlabel("Distance from divide (km)", fontsize=fs)
    fig.savefig("../plots/figure_11.pdf", dpi=300)


def fabric(fn, xlim=(-30, 30), ylim=(-250, 3250), clim=(-0.25, 0.25)):
    vtu = fastvtulib.get_structured_vtu(fn)
    vtu.rawdata_dict["xvel"] = vtu.rawdata_dict["aiflow"][:, 0]
    check_offset_guesses = np.arange(-10000, 10000, 10.0)
    xvels = vtu.get_pts_2d(["xvel"], check_offset_guesses, np.ones_like(check_offset_guesses) * ylim[1] * 0.9)["xvel"]
    offset = check_offset_guesses[np.argmin(np.abs(xvels))] / 1000.0

    tri = Triangulation(np.array([rc[0] / 1000.0 for rc in vtu.raw_coords[:, 0]]), np.array([rc[1] for rc in vtu.raw_coords[:, 0]]), vtu.simptt)

    fig, ((ax1, cax1), (ax2, cax2), (ax3, cax3)) = plt.subplots(figsize=(7.047, 5.0), nrows=3, ncols=2, gridspec_kw={"left": 0.1, "right": 0.92, "bottom": 0.09, "top": 0.95, "width_ratios": (1, 0.02), "wspace": 0.05, "hspace": 0.06})
    for ax in (ax1, ax2, ax3):
        ax.set_facecolor("k")

    cm1 = ax1.tricontourf(tri, vtu.rawdata_dict["aiflow"][:, 0], cmap="RdBu_r", levels=np.linspace(clim[0], clim[1], 256), extend="both")
    cm2 = ax2.tricontourf(tri, vtu.rawdata_dict["tensorfabric 2"], cmap="BrBG", norm=anisolib.fab_norm, levels=np.linspace(0, 1, 101), extend="neither")
    rot = fabricplotlib.fabric_to_ver_rot(vtu.rawdata_dict["tensorfabric 1"], vtu.rawdata_dict["tensorfabric 2"], vtu.rawdata_dict["tensorfabric 3"])
    cm3 = ax3.tricontourf(tri, rot, levels=np.linspace(-45, 45, 256), cmap="PuOr", extend="both")

    levels = [5000, 10000, 25000, 50000, 100000, 250000]
    CS = ax1.tricontour(tri, vtu.rawdata_dict["age"], colors="k", levels=levels, linewidths=0.5)
    fmt = {level: "{:d} ka".format(level // 1000) for level in levels}
    ax1.clabel(CS, CS.levels, inline=True, fmt=fmt, fontsize=8, manual=[(28, 2900), (24, 2800), (24, 2500), (24, 2000), (24, 3200), (24, 500)])
    for ctr in CS.collections:
        for path in ctr.get_paths():
            v = path.vertices
            x = v[:, 0]
            y = v[:, 1]
            mask = np.abs(x) < 5
            if np.any(mask):
                print("Bump height:", np.max(y[mask]) - np.min(y[mask]))

    for c in cm1.collections:
        c.set_edgecolor("face")

    cbr1 = plt.colorbar(cm1, cax=cax1, orientation="vertical", ticks=(clim[0], 0, clim[1]))
    cbr1.set_label(label=r"$v_x$ (m yr$^{-1})$", size=fs)
    cbr1.ax.tick_params(axis="both", which="major", labelsize=fs)

    cbr2 = plt.colorbar(cm2, cax=cax2, orientation="vertical", ticks=(0, 1 / 3, 2 / 3, 1.0))
    cbr2.set_label(label=r"$a^{(2)}_{zz}$", size=fs)
    cbr2.ax.set_yticklabels(["0", r"$\frac{1}{3}$", r"$\frac{2}{3}$", "1"])

    cbr3 = plt.colorbar(cm3, cax=cax3, orientation="vertical", ticks=(-45, 0, 45))
    cbr3.set_label(label="Rotation from vertical (\N{DEGREE SIGN})", size=fs)

    for letter, ax in zip("abcdefghijkl", (ax1, ax2, ax3)):
        ax.text(0.01, 0.82, letter, transform=ax.transAxes, fontsize=bfs)
        ax.tick_params(axis="both", which="major", labelsize=fs)
        ax.set_ylim(*ylim)
        # ax.set_xlim(*xlim)
        ax.set_ylabel("Elevation (m)", fontsize=fs)
        if abs(offset) > 0.1:
            if False:  # ax is ax3:
                label = "Modeled divide"
            else:
                label = None
            ax.axvline(offset, linestyle="dotted", color="k", linewidth=0.75, label=label)

    if abs(offset) > 0.1:
        # trans = transforms.blended_transform_factory(ax3.transData, ax3.transAxes)
        # ax3.text(offset, 0.5, 'Modeled divide', rotation=90, ha='right', va='center', transform=trans, fontsize=fs)
        labelLines(ax3.get_lines(), zorder=9, backgroundcolor="none", fontsize=fs, va="top", yoffsets=3000)
        # labelLines(ax2.get_lines(), zorder=9, fontsize=fs, yoffsets=[2500] + [3000 for y in site_data['use_sites'][1:]])

    for ax in (ax1, ax2):
        ax.tick_params(axis="x", which="both", labelbottom=False, labelleft=False)

    return fig, (ax1, ax2, ax3)


def main():
    divide_2p("")
    domec_fabric()


if __name__ == "__main__":
    main()
