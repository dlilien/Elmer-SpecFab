#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.


"""

"""

import glob

import numpy as np
import matplotlib.pyplot as plt
from matplotlib import gridspec

from modeltools.lib import fastvtulib

FS = 10
BFS = 14


def edc_plot_calib(fabname="tensorfabric", tube="", subsr="", folder="domec"):
    y_pts = 500
    edc_fn = "../edc_data/EDC_fabric_durand.txt"
    edc_ds = np.genfromtxt(edc_fn, names=True)

    gs = gridspec.GridSpec(1, 1, left=0.19, right=0.98, top=0.98, bottom=0.13)
    fig = plt.figure(figsize=(3.4, 4.0))
    ax = fig.add_subplot(gs[0, 0])

    cm = plt.get_cmap("tab20")
    ax.plot(edc_ds["a1"], edc_ds["depth"], marker="s", markerfacecolor=cm(1), markeredgecolor="none", linestyle="none", markersize=5)
    ax.plot(edc_ds["a2"], edc_ds["depth"], marker="s", markerfacecolor=cm(3), markeredgecolor="none", linestyle="none", markersize=5)
    ax.plot(edc_ds["a3"], edc_ds["depth"], marker="s", markerfacecolor=cm(5), markeredgecolor="none", linestyle="none", markersize=5)

    ax.set_xlim(0, 1)
    ax.set_xticks((0, 1.0 / 3.0, 2.0 / 3.0, 1))
    ax.set_xticklabels(["0", r"$\frac{1}{3}$", r"$\frac{2}{3}$", "1"], fontsize=FS)
    H = 3250
    ax.set_ylim(H, 0)
    ax.set_yticks(np.flipud(np.arange(0, H, 500)))
    ax.set_yticks(np.flipud(np.arange(0, H, 250)), minor=True)
    ax.set_ylabel("Depth (m)", fontsize=FS)
    ax.set_xlabel(r"$\lambda_i$", fontsize=FS)
    ax.tick_params(axis="both", which="major", labelsize=FS)
    ls = [["solid"], ["dashed"]]
    ns = [["Lab"], ["Ice core"]]
    for calib, names, linestyles in zip(["", "_calib"], ns, ls):
        fn_glob = "../domec_transect/domec_L6/domec{:s}_full{:s}{:s}_t????.vtu".format(tube, subsr, calib)

        files = glob.glob(fn_glob)
        inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))

        def plt_frame(frame, set_time=True):
            fn = files[frame]
            vtu = fastvtulib.get_structured_vtu(fn)
            dists = [0.0]

            for dist, linestyle, name in zip(dists, linestyles, names):
                x = np.ones((y_pts,)) * dist
                y = np.linspace(np.min(vtu.coords[:, 1]), np.max(vtu.coords[:, 1]), y_pts)
                test_p = vtu.get_pts_2d(["{:s} 1".format(fabname)], x, y)
                y_valid = y[~np.isnan(test_p["{:s} 1".format(fabname)])]
                y = np.linspace(np.min(y_valid), np.max(y_valid), y_pts)[::-1]
                depths = y[0] - y
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

                evals, evecs = np.linalg.eig(A)
                evals = np.sort(evals)
                ax.plot(evals[:, 2], depths, color="C0", linestyle=linestyle, zorder=99)
                ax.plot(evals[:, 1], depths, color="C1", linestyle=linestyle, zorder=99)
                ax.plot(evals[:, 0], depths, color="C2", linestyle=linestyle, zorder=99)

        plt_frame(inds[-1], set_time=False)

    ax.plot([], [], linestyle="none", marker="s", markerfacecolor="k", markeredgecolor="none", label="EDC data")
    for linestyle, name in zip(ls, ns):
        ax.plot([], [], linestyle=linestyle[0], color="k", label=name[0])
    ax.legend(loc="upper right", fontsize=FS, frameon=False)

    fig.savefig("../plots/figure_6.pdf", dpi=300)


def main():
    edc_plot_calib(subsr="_aniso", tube="_tube")


if __name__ == "__main__":
    main()
