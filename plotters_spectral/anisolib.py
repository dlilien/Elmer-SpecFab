# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © D. A. Lilien <david.lilien@umanitoba.ca>, 2022
#
# Distributed under terms of the GNU GPL3.0 license.


"""

"""

import numpy as np
from matplotlib import colors

fab_norm = colors.TwoSlopeNorm(vmin=0.0, vcenter=0.333, vmax=1)
cartoon_norm = colors.LogNorm(vmin=0.01, vmax=1.0, clip=False)
fabs = ["fabric 1", "fabric 2", "fabric 3", "fabric 4", "fabric 5"]


def fabric_dict_to_a2(fabric_dict, twod=True):
    return elmer_fabric_to_a2(*[fabric_dict[fab] for fab in fabs], twod=twod)


def elmer_fabric_to_a2(a11, a22, a12, a23, a13, twod=True):
    a33 = 1.0 - a11 - a22
    a33[np.isnan(a11)] = np.nan
    a33[a33 < 0.0] = 0.0
    if twod:
        return np.array([[a11, a13, a12], [a13, a33, a23], [a13, a23, a22]])
    else:
        return np.array([[a11, a12, a13], [a12, a22, a23], [a13, a23, a33]])


def elmer_fabric_to_eigs(a11, a22, a12, a23, a13, prec=1.0e-8):
    """Somehow this is introducing NaNs in places where we should have values...not yet sure why"""
    A = elmer_fabric_to_a2(a11, a22, a12, a23, a13)
    if len(A.shape) == 3:
        A = np.transpose(A, (2, 0, 1))
    A[np.abs(A) < prec] = 0.0
    mask = np.all(np.all(~np.isnan(A), axis=-1), axis=-1)
    print("Starting eigs")
    w, v = np.linalg.eigh(A[mask, :, :])
    print("eigs found")
    w = np.abs(w)

    # force all of the eigenvectors to be in the -z direction
    if len(A.shape) == 2:
        for i in range(len(w)):
            if v[2, i] < 0:
                v[:, i] = -v[:, i]
    else:
        for j in range(w.shape[0]):
            for i in range(w.shape[1]):
                if v[j, 2, i] < 0:
                    v[j, :, i] = -v[j, :, i]
    if not np.all(mask):
        v_out = np.ones_like(A)
        w_out = np.ones(A.shape[:-1])
        if len(A.shape) == 2:
            pass
        else:
            v_out[:, :, :] = np.nan
            v_out[mask, :, :] = v
            w_out[:, :] = np.nan
            w_out[mask, :] = w
    else:
        w_out = w
        v_out = v
    return w_out, v_out


def v_to_angles(v):
    sigma_x = np.arctan2(v[:, 2, 1], v[:, 2, 2]) * 180.0 / np.pi
    sigma_y = np.arctan2(-v[:, 2, 0], np.sqrt(v[:, 2, 1] ** 2.0 + v[:, 2, 2] ** 2.0)) * 180.0 / np.pi
    sigma_z = np.arctan2(v[:, 0, 1], v[:, 0, 0]) * 180.0 / np.pi
    return sigma_x, sigma_y, sigma_z


def naive_vectors(w, v, num):
    weight1 = np.atleast_2d(np.random.normal(size=num) / w[0] ** 2.0).transpose()
    weight2 = np.atleast_2d(np.random.normal(size=num) / w[1] ** 2.0).transpose()
    weight3 = np.atleast_2d(np.random.normal(size=num) / w[2] ** 2.0).transpose()

    return weight1 * v[:, 0] + weight2 * v[:, 1] + weight3 * v[:, 2]


def woodcock_param(w, tol=1.0e-2, maxk=100.0):
    k = np.log10(np.abs(w[:, 2] / w[:, 1])) / np.log10(np.abs(w[:, 1] / w[:, 0]))
    # print(w[:, 2], w[:, 1], w[:, 0])
    # k[np.abs(1. - w[:, 1] / w[:, 0]) < tol] = tol
    # k[np.abs(1. - w[:, 2] / w[:, 1]) < tol] = 0
    # if they are all the same, go with 1
    # k[np.logical_and(np.abs(1. - w[:, 2] / w[:, 1]) < tol, np.abs(1. - w[:, 1] / w[:, 0]) < tol)] = 1
    k[np.logical_or(np.isnan(k), np.isinf(k))] = maxk
    k[np.logical_and(np.abs(w[:, 1]) < tol, np.abs(w[:, 0]) < tol)] = maxk
    return k
