# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © N. M. Rathmann <rathmann@nbi.ku.dk> and D. A. Lilien <david.lilien@umanitoba.ca>, 2022
#
# Distributed under terms of the GNU GPL3.0 license.


import os
import pickle
import numpy as np
import quaternion as qt  # pip install numpy-quatern
import pandas as pd

from specfabpy import specfabpy as sf

from header import get_v_angles

L = 20  # Spectral truncation
pickle_fn = "experiments.p"
RELOAD = False


class Experiment:
    def __init__(self, name, src, fn, color, temp, tot_strain, strain_type="uc", skip_rows=0, sep=","):
        self.name = name
        if src not in ["Fan", "Qi", "Hunter"]:
            raise ValueError("Unsupported source")
        self.src = src
        self.color = color
        self.temp = temp
        self.tot_strain = tot_strain
        if strain_type not in ["uc", "cc", "ut", "ss"]:
            raise ValueError("Unsupported strain type")
        self.strain_type = strain_type

        self.load(fn, skip_rows, sep)

    def load(self, fn, skip_rows, sep):
        df = pd.read_csv(fn, skiprows=skip_rows, sep=sep)

        # For SPICE files remove rows with "unrecognized grains"
        m = ~df.astype(str).apply(lambda x: x.str.contains("grain", regex=False)).any(axis=1)
        df = df[m]

        # Processed EBSD data in MTEX (MATLAB) is saved as quaternion coords.
        qs_comps = df.to_numpy()[:, :]
        qs = qt.as_quat_array(qs_comps)
        sphcoords = qt.as_spherical_coords(qs)
        # COLAT [0;pi], LON [0;2pi]
        qcolat, qlon = sphcoords[:, 0], sphcoords[:, 1]
        # qcolat, qlon = np.array([np.pi/2 * 1,]), np.array([np.pi/2 * 1,]) # debug

        # Add reflected vectors; ODFs are antipodally symmetric, so doesn't change statistics
        if 1:
            qcolat = np.hstack((qcolat, np.pi - qcolat))
            qlon = np.hstack((qlon, qlon - np.pi))

        # Determine ODF from a2 and a4

        lm, nlm_len = sf.init(4)
        # The expansion coefficients
        nlm_L4 = np.zeros((nlm_len), dtype=np.complex64)
        nlm_L2 = nlm_L4.copy()
        caxes = np.array([[np.cos(p) * np.sin(t), np.sin(p) * np.sin(t), np.cos(t)] for t, p in zip(qcolat, qlon)])
        a2 = np.array([np.einsum("i,j", c, c) for c in caxes]).mean(axis=0)
        a4 = np.array([np.einsum("i,j,k,l", c, c, c, c) for c in caxes]).mean(axis=0)
        #    print(a2)
        nlm_L2[:6] = sf.a2_to_nlm(a2)
        nlm_L4[:16] = sf.a4_to_nlm(a4)
        #    nlm_L4[6:] = 0.0

        # Rotated frame

        # sym. axis = largest eig dir. for single max, smallest eig dir. for girdle
        Ilami = 0 if self.strain_type == "ut" else 2
        (v1_colat, v1_lon, _) = get_v_angles(a2, Ilami=Ilami)
        # v1_colat = 0 # only horiz rotation (for debugging)

        print("true v1 colat, lon = %f, %f (deg.) " % (np.rad2deg(v1_colat), np.rad2deg(v1_lon)))
        # THIS is the nlm array from which spectral coefs are derived for correlation
        nlmr_L4 = sf.rotate_nlm4(nlm_L4, 0, -v1_lon)
        nlmr_L4 = sf.rotate_nlm4(nlmr_L4, +v1_colat, 0)
        self.a2 = sf.a2(nlmr_L4)
        self.a4 = sf.a4(nlmr_L4)


def load_pickle(RELOAD):
    if not os.path.exists(pickle_fn) or RELOAD:
        PIL135 = Experiment("PIL135", "Qi", "latrot-validation/Qi_etal/PIL135 shear 20 clean.ctf.csv", "tab:olive", -30.5, 2.6, strain_type="ss", skip_rows=0, sep=",")

        PIL94 = Experiment("PIL94", "Qi", "latrot-validation/Qi_etal/PIL94 shear 7 clean.ctf.csv", "tab:olive", -5.2, 1.5, strain_type="ss", skip_rows=0, sep=",")

        PIL255 = Experiment("PIL268", "Qi", "latrot-validation/Fan_etal/series_-30/PIL268_30mu.ctf.csv", "tab:red", -20.0, 0.20, strain_type="uc", skip_rows=0, sep=",")

        PIL007 = Experiment("PIL007", "Qi", "latrot-validation/Fan_etal/series_-10/PIL007_30mu.crc.csv", "tab:red", -10.0, 0.19, strain_type="uc", skip_rows=0, sep=",")

        PIL268 = Experiment("PIL255", "Fan", "latrot-validation/Fan_etal/series_-20/PIL255_30mu.ctf.csv", "tab:red", -30.0, 0.21, strain_type="uc", skip_rows=0, sep=",")

        D53 = Experiment("D5-3", "Hunter", "latrot-validation/Hunter_etal/D5-3/D5-3_pf.csv", "tab:gray", -1.0, 0.4, strain_type="uc", skip_rows=0, sep=",")

        D51 = Experiment("D5-1", "Hunter", "latrot-validation/Hunter_etal/D5-1/D5-1_pf.csv", "tab:gray", -7.0, 0.4, strain_type="uc", skip_rows=0, sep=",")

        MD22 = Experiment("MD22", "Hunter", "latrot-validation/Hunter_etal/MD22/md22.csv", "tab:gray", -7.0, 0.4, strain_type="uc", skip_rows=0, sep=",")

        experiments = [PIL135, PIL94, PIL268, PIL255, PIL007, MD22, D53, D51]
        with open(pickle_fn, "wb") as fout:
            pickle.dump(experiments, fout)
    else:
        with open(pickle_fn, "rb") as fin:
            experiments = pickle.load(fin)
    return experiments


if __name__ == "__main__":
    load_pickle(RELOAD)
