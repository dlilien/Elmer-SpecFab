#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.


"""
Make Figure 2 and supplementary figure 2
"""
import h5py
import glob
import yaml

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec, ticker

from modeltools.lib import fastvtulib

π = np.pi
FS = 10
BFS = 15

edc_fn = "../edc_data/EDC_fabric_durand.txt"
edc_ds = np.genfromtxt(edc_fn, names=True)


def modeled_vs_measured_pres():
    # We want to find the filename for the highest index where we have all 3
    with open("pres_sites.yml") as fin:
        site_data = yaml.safe_load(fin)

    # From paper: The angle α is measured by compass with ±15° uncertainty for georeferencing
    # the data. Here, we use polar stereographic coordinates where anticlockwise
    # rotation is positive.

    fn_glob = "../domec_transect/domec_v2/domec_tube_full_aniso_calib_t????.vtu"
    files = glob.glob(fn_glob)
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    fn = files[inds[-1]]
    calib_vtu = fastvtulib.get_structured_vtu(fn)

    fn_glob = "../domec_transect/domec_v2/domec_tube_full_aniso_t????.vtu"
    files = glob.glob(fn_glob)
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    fn = files[inds[-1]]
    lab_vtu = fastvtulib.get_structured_vtu(fn)
    y_pts = 1000

    figdum, axdum = plt.subplots()
    ncol = 7
    gs = gridspec.GridSpec(nrows=1, ncols=ncol, left=0.09, top=0.97, bottom=0.24, right=0.995, wspace=0.15, width_ratios=(1, 1, 1, 1, 1, 1, 1))
    fig = plt.figure(figsize=(7.8, 3.2))
    top_axes = [fig.add_subplot(gs[0, i]) for i in range(ncol)]
    top_axes[3].plot(edc_ds["a2"] - edc_ds["a3"], edc_ds["depth"], marker="s", color="darkred", linestyle="none", markersize=4, label="Ice-core data")
    top_axes[0].plot([], [], marker="s", color="darkred", linestyle="none", markersize=4, label="Ice-core data")
    for site, ax_da in zip(site_data["sparse_sites"], top_axes):
        measured_data = h5py.File("../edc_data/ApRES_DomeC/{:s}_ha.mat".format(site[5:].replace(".", "d")))
        # Reza uses a "true north" coordinate system where he has flow at 45
        # - float(angle_data['sites'][site]['azimuth'])
        # print(float(angle_data['sites'][site]['azimuth']))
        ax_da.plot(measured_data["ha"][:].flatten(), measured_data["Z"][:].flatten(), linestyle="None", marker="o", mfc="none", mec="k", markersize=4, label="pRES")

        x = float(site_data["sites"][site]["dist"])
        plot_model(lab_vtu, ax_da, axdum, x, y_pts, linestyle="dashed", color="0.6", label="Lab. cal. model")
        plot_model(calib_vtu, ax_da, axdum, x, y_pts, color="0.2", linestyle="solid", label="Ice-core cal. model")

        ax_da.set_xlim(0, 0.15)

    for letter, tax, name in zip("abcdefghijklmnop", top_axes, site_data["sparse_sites"] + site_data["sparse_sites"]):
        tax.set_ylim(2000, 0)
        tax.tick_params(axis="both", which="major", labelsize=FS)

        # tax.text(0.05, 0.90, '{:s}1'.format(letter), transform=tax.transAxes, fontsize=BFS)
        # tax.text(0.35, 0.90, name[5:], transform=tax.transAxes, fontsize=FS)
        # bax.text(0.05, 0.90, '{:s}2'.format(letter), transform=bax.transAxes, fontsize=BFS)
        # bax.text(0.35, 0.90, name[5:], transform=bax.transAxes, fontsize=FS)

        tax.text(0.5, 0.92, "{:s}".format(letter), ha="center", transform=tax.transAxes, fontsize=BFS)
        # tax.text(0.35, 0.90, name[5:], transform=tax.transAxes, fontsize=FS)
        # bax.text(0.35, 0.90, name[5:], transform=bax.transAxes, fontsize=FS)
        # tax.set_title(name[5:], fontsize=12)

    for ax in top_axes[1:]:
        ax.yaxis.set_major_formatter(ticker.NullFormatter())
    top_axes[3].set_xlabel("Horizontal eigenvalue difference", fontsize=10)
    top_axes[0].set_ylabel("Depth (m)", fontsize=10)
    top_axes[3].legend(loc="upper left", ncol=4, bbox_to_anchor=(-3.0, -0.18), fontsize=10)

    fig.savefig("../plots/figure_7.pdf")


def plot_model(vtu, ax_dl, ax_theta, x, y_pts, fabname="tensorfabric", linestyle="solid", color="k", label=None):
    x = np.ones((y_pts,)) * x
    y = np.linspace(np.min(vtu.coords[:, 1]), np.max(vtu.coords[:, 1]), y_pts)
    test_p = vtu.get_pts_2d(["{:s} 1".format(fabname)], x, y)
    y_valid = y[~np.isnan(test_p["{:s} 1".format(fabname)])]
    y = np.linspace(np.min(y_valid), np.max(y_valid), y_pts)[::-1]
    pts = vtu.get_pts_2d(["{:s} {:d}".format(fabname, i + 1) for i in range(5)] + ["temp", "temp homologous"], x, y)
    A = np.zeros((x.shape[0], 3, 3))
    A[:, 0, 0] = pts["{:s} 1".format(fabname)]
    A[:, 1, 1] = 1.0 - pts["{:s} 1".format(fabname)] - pts["{:s} 2".format(fabname)]
    A[:, 2, 2] = pts["{:s} 2".format(fabname)]
    A[:, 0, 1] = pts["{:s} 5".format(fabname)]
    A[:, 1, 0] = pts["{:s} 5".format(fabname)]
    A[:, 2, 1] = pts["{:s} 4".format(fabname)]
    A[:, 1, 2] = pts["{:s} 4".format(fabname)]
    A[:, 2, 0] = pts["{:s} 3".format(fabname)]
    A[:, 0, 2] = pts["{:s} 3".format(fabname)]

    w, v = np.linalg.eig(A)  # [:, :2, :2])
    order = np.argsort(w, axis=1)
    w = np.take_along_axis(w, order, 1)
    # v = np.take_along_axis(v, order, 2)
    da = w[:, 1] - w[:, 0]
    theta = np.arctan2(v[:, 0, 0], v[:, 1, 0]) * 180.0 / np.pi
    theta[theta < 0.0] = theta + 180.0
    depth = np.max(y) - y
    ax_dl.plot(A[:, 1, 1] - A[:, 0, 0], depth, color=color, linestyle=linestyle, label=label)
    ax_theta.plot(theta, y, color=color, linestyle=linestyle, label=label)


def modeled_vs_measured_pres_supp():
    # We want to find the filename for the highest index where we have all 3
    with open("pres_sites.yml") as fin:
        site_data = yaml.safe_load(fin)

    # From paper: The angle α is measured by compass with ±15° uncertainty for georeferencing
    # the data. Here, we use polar stereographic coordinates where anticlockwise
    # rotation is positive.

    fn_glob = "../domec_transect/domec_v2/domec_tube_full_aniso_calib_t????.vtu"
    files = glob.glob(fn_glob)
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    fn = files[inds[-1]]
    calib_vtu = fastvtulib.get_structured_vtu(fn)

    fn_glob = "../domec_transect/domec_v2/domec_tube_full_aniso_t????.vtu"
    files = glob.glob(fn_glob)
    inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))
    fn = files[inds[-1]]
    lab_vtu = fastvtulib.get_structured_vtu(fn)
    y_pts = 1000

    all_sites = [key for key in site_data["sites"].keys()]
    for key in all_sites[::-1]:
        if key[0] != "p":
            all_sites.remove(key)
        if key.lower() in [k.lower() for k in site_data["sparse_sites"]]:
            all_sites.remove(key)
        if key.lower() == "pres_w3.5":
            all_sites.remove(key)
        if key.lower() == "pres_w9.0":
            all_sites.remove(key)

    nums = np.array([float(site.split("_")[1].replace("E", "-").replace("W", "")) for site in all_sites])
    order = np.argsort(nums)

    figdum, axdum = plt.subplots()
    ncol = 6
    gs = gridspec.GridSpec(nrows=2, ncols=ncol, left=0.1, top=0.99, bottom=0.17, right=0.995, wspace=0.15, hspace=0.05)
    fig = plt.figure(figsize=(7.047, 5.75))
    top_axes = [fig.add_subplot(gs[0, i]) for i in range(ncol)]
    bot_axes = [fig.add_subplot(gs[1, i]) for i in range(ncol)]
    for i, ax_da in zip(order, top_axes + bot_axes):
        site = all_sites[int(i)]
        measured_data = h5py.File("../edc_data/ApRES_DomeC/{:s}_ha.mat".format(site[5:].replace(".", "d")))
        # Reza uses a "true north" coordinate system where he has flow at 45
        # theta = measured_data['v2'][:].flatten() - 45.0
        # - float(angle_data['sites'][site]['azimuth'])
        # print(float(angle_data['sites'][site]['azimuth']))
        ax_da.plot(measured_data["ha"][:].flatten(), measured_data["Z"][:].flatten(), linestyle="None", marker="o", mfc="none", mec="k", markersize=4, label="pRES")

        x = float(site_data["sites"][site]["dist"])
        plot_model(lab_vtu, ax_da, axdum, x, y_pts, linestyle="dashed", color="0.6", label="Lab. cal. model")
        plot_model(calib_vtu, ax_da, axdum, x, y_pts, color="0.2", linestyle="solid", label="Ice-core cal. model")

        ax_da.set_xlim(0, 0.15)
        ax_da.set_xticks([0, 0.1])

    for letter, tax, sitenum in zip("abcdefghijklmnop", top_axes + bot_axes, order):
        site = all_sites[int(sitenum)]
        tax.set_ylim(2000, 0)
        tax.tick_params(axis="both", which="major", labelsize=FS)

        # tax.text(0.05, 0.90, '{:s}1'.format(letter), transform=tax.transAxes, fontsize=BFS)
        # tax.text(0.35, 0.90, name[5:], transform=tax.transAxes, fontsize=FS)
        # bax.text(0.05, 0.90, '{:s}2'.format(letter), transform=bax.transAxes, fontsize=BFS)
        # bax.text(0.35, 0.90, name[5:], transform=bax.transAxes, fontsize=FS)

        tax.text(0.5, 0.92, "{:s}".format(letter), ha="right", transform=tax.transAxes, fontsize=BFS)
        tax.text(0.5, 0.92, " " + site[5:], ha="left", transform=tax.transAxes, fontsize=FS)
        # tax.text(0.35, 0.90, name[5:], transform=tax.transAxes, fontsize=FS)
        # bax.text(0.35, 0.90, name[5:], transform=bax.transAxes, fontsize=FS)
        # tax.set_title(name[5:], fontsize=12)

    for ax in top_axes[1:] + bot_axes[1:]:
        ax.yaxis.set_major_formatter(ticker.NullFormatter())
    for ax in top_axes:
        ax.xaxis.set_major_formatter(ticker.NullFormatter())
    bot_axes[3].set_xlabel("Horizontal eigenvalue difference                                   ", fontsize=10)
    top_axes[0].set_ylabel("Depth (m)", fontsize=10)
    bot_axes[0].set_ylabel("Depth (m)", fontsize=10)
    bot_axes[0].legend(loc="upper left", ncol=4, bbox_to_anchor=(0.0, -0.18), fontsize=8)

    fig.savefig("../plots/supplementary_figure_2.pdf")


def main():
    modeled_vs_measured_pres()
    modeled_vs_measured_pres_supp()


if __name__ == "__main__":
    main()
