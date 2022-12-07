#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Make plots of time-dependent quantities at EDC (Figure 4 of the paper).
"""
import glob

import numpy as np
import matplotlib.pyplot as plt
import matplotlib

matplotlib.rcParams["contour.negative_linestyle"] = "solid"

from modeltools.lib import fastvtulib

FS = 8
BFS = 12


def main(time_offset=250000):
    names = {"Ice core": "aniso_calib", "Lab": "aniso"}
    offset_dict = {}
    time_dict = {}
    height_dict = {}
    for name, fn_name in names.items():
        fn_glob = "../domec_transect/domec_L6/domec_tube_full_{:s}_t????.vtu".format(fn_name)
        files = glob.glob(fn_glob)
        inds = np.argsort(np.array([float(fn[-8:-4]) for fn in files]))

        offset_times = time_offset - np.hstack(([0.0, 1.25], np.arange(1, len(files) - 1) * 100.0))
        offsets = np.zeros_like(offset_times)
        heights = np.zeros_like(offset_times)

        for i, ind in enumerate(inds):
            fn = files[ind]
            vtu = fastvtulib.get_structured_vtu(fn)
            vtu.rawdata_dict["xvel"] = vtu.rawdata_dict["aiflow"][:, 0]
            check_offset_guesses = np.arange(-10000, 10000, 10.0)
            xvels = vtu.get_pts_2d(["xvel"], check_offset_guesses, np.ones_like(check_offset_guesses) * 3000.0)["xvel"]
            offsets[i] = check_offset_guesses[np.argmin(np.abs(xvels))] / 1000.0
            heights[i] = np.max(vtu.coords[:, 1]) - np.min(vtu.coords[vtu.coords[:, 0] < 100.0, 1])
        offset_dict[name] = offsets
        height_dict[name] = heights
        time_dict[name] = offset_times

    acc_arr = np.genfromtxt("../edc_data/AICC2012_acc.csv", skip_header=1)
    acc_times = acc_arr[:, 0]
    acc = acc_arr[:, 1]
    temp_arr = np.genfromtxt("../edc_data/jouzel_2007_EDC_dD_temp.tab", skip_header=1)
    temp_times = temp_arr[:, 1] * 1000.0
    temp = temp_arr[:, 3] - 54.0

    fig, axes = plt.subplots(4, 1, sharex=True, figsize=(5.5, 3.25))
    (temp_ax, acc_ax, divide_ax, thick_ax) = axes
    acc_ax.plot(acc_times / 1000.0, acc * 100, color="k")
    acc_ax.set_ylabel("Acc. rate\n(cm yr$^{-1}$)", fontsize=FS)
    acc_ax.set_ylim(1, 5)
    acc_ax.set_yticks([1, 2, 3, 4, 5])

    temp_ax.plot(temp_times / 1000.0, temp, color="k")
    temp_ax.set_ylabel("Temp.\n(\N{DEGREE SIGN}C)", fontsize=FS)
    temp_ax.set_ylim(-65, -45)
    temp_ax.set_yticks([-65, -60, -55, -50, -45])

    for name, ls in zip(names, ["solid", "dashed", "dotted"]):
        divide_ax.plot(time_dict[name] / 1000.0, offset_dict[name], color="0.6", ls=ls, label=name)
        thick_ax.plot(time_dict[name] / 1000.0, height_dict[name] / 1000.0, color="0.6", ls=ls, label=name)

    divide_ax.set_ylabel("Divide\nposition\n(km)", fontsize=FS)
    divide_ax.set_ylim(-1, 1)
    divide_ax.set_yticks([-1, 0, 1])

    divide_ax.legend(loc="upper right", frameon=False, fontsize=FS)

    thick_ax.set_ylabel("Ice thick.\n(km)", fontsize=FS)
    thick_ax.set_ylim(3.53, 3.55)
    thick_ax.set_yticks([3.53, 3.54, 3.55])

    axes[-1].set_xlabel("Time (kyr before 1950)", fontsize=FS)

    for letter, ax in zip("abcdefghijkl", axes):
        ax.text(0.01, 0.99, letter, ha="left", va="top", transform=ax.transAxes, fontsize=BFS)
        ax.tick_params(axis="both", which="major", labelsize=FS)
        ax.set_xlim(250, 0)

    fig.tight_layout(pad=0.3)
    fig.savefig("../plots/figure_4.pdf", dpi=300)


if __name__ == "__main__":
    main()
