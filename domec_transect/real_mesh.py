#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.
import sys

import subprocess as sp
import numpy as np
from scipy.interpolate import LinearNDInterpolator
import shapely.geometry

import xarray as xr
import rioxarray
from scipy.ndimage import gaussian_filter, gaussian_filter1d
import matplotlib.pyplot as plt

import osr
from modeltools.lib.geo_classes import Interpolator, Centerline
from modeltools.lib import glib
from modeltools.lib import gmshlib

from pygeotools.lib.geolib import ll2sps

rioxarray
V_SCALE = 4.0
RES = 300.0

debug = False

data_dir = sys.argv[1]

EDC = (1359695.36, -894852.43)
x_edc = EDC[0]
y_edc = EDC[1]
FIRN_DIP = 18.0
BETA_TAPER = 2000.0
basal_melt_rate = 0.01


def get_conversion(t_srs):
    out_cs = osr.SpatialReference()
    out_cs.SetFromUserInput(t_srs)
    try:
        out_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except AttributeError:
        pass

    wgs84_cs = out_cs.CloneGeogCS()
    wgs84_cs.ExportToPrettyWkt()
    try:
        wgs84_cs.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    except AttributeError:
        pass

    transform_WGS84_To_srs = osr.CoordinateTransformation(out_cs, wgs84_cs)
    return transform_WGS84_To_srs.TransformPoints, out_cs.ExportToPrettyWkt()


sps2ll = get_conversion("EPSG:3031")[0]


def interpolate(x, y, x0, y0, q):
    """
    Interpolate bilinearly.

    Parameters
    ----------
    x : np.ndarray
        x coordinates for data matrix q.
    y : TYPE
        y coordinates for data matrix q.
    x0 : TYPE
        The x value of point at which to interpolate.
    y0 : TYPE
        The y value of the point at which to interpolate.
    q : TYPE
        The data matrix.

    Returns
    -------
    p : float
        The interpolated value.

    """
    dx = x[1] - x[0]
    dy = y[1] - y[0]

    i = max(min(int((y0 - y[0]) / dy), len(y) - 1), 0)
    j = max(min(int((x0 - x[0]) / dx), len(x) - 1), 0)

    p = q[i, j] + (x0 - x[j]) * (q[i, j + 1] - q[i, j]) / dx + (y0 - y[i]) * (q[i + 1, j] - q[i, j]) / dy + (x0 - x[j]) * (y0 - y[i]) * (q[i + 1, j + 1] + q[i, j] - q[i + 1, j] - q[i, j + 1]) / (dx * dy)
    return p


def streamline(x, y, vx, vy, x0, y0, stepsize=50.0, minvel=0.01, timesize=None, maxpts=1e4, trackbetter=True):
    """x,y are coords of vx and vy, x0,y0 is starting point."""
    u = interpolate(x, y, x0, y0, vx)
    v = interpolate(x, y, x0, y0, vy)

    speed = np.sqrt(u**2 + v**2)

    X = [x0]
    Y = [y0]
    U = [speed]
    dts = [0.0]

    k = 0

    while speed > minvel:
        k = k + 1
        if timesize is None:
            dt = stepsize / speed
        else:
            dt = timesize

        # To track better, take lots of little steps to form the flowline
        if trackbetter:
            extra_steps = 100
            dt = dt / float(extra_steps)
            for i in range(extra_steps):
                x0 = x0 + dt * u
                y0 = y0 + dt * v
                u = interpolate(x, y, x0, y0, vx)
                v = interpolate(x, y, x0, y0, vy)
                if np.isnan(u):
                    break
            if np.isnan(u):
                break
        else:
            x0 = x0 + dt * u
            y0 = y0 + dt * v

        X.append(x0)
        Y.append(y0)
        U.append(speed)
        dts.append(dt)

        u = interpolate(x, y, x0, y0, vx)
        v = interpolate(x, y, x0, y0, vy)
        speed = np.sqrt(u**2 + v**2)
        if (speed > 1.0e5) or np.isnan(speed):
            if k < 100:
                print("eek")
            break
        if k > maxpts:
            print("Iters capped")
            break

    t = np.cumsum(dts)
    return X, Y, U, t


def make_centerlines(x, y, vx, vy, xs, ys, maxpts_down=750, maxpts_up=45, mv=0.01):
    """
    Go down and back from the starting point.

    Parameters
    ----------
    x : TYPE
        DESCRIPTION.
    y : TYPE
        DESCRIPTION.
    vx : TYPE
        DESCRIPTION.
    vy : TYPE
        DESCRIPTION.
    xs : TYPE
        DESCRIPTION.
    ys : TYPE
        DESCRIPTION.
    maxpts : TYPE, optional
        DESCRIPTION. The default is 1e4.
    mv : TYPE, optional
        DESCRIPTION. The default is 0.1.

    Returns
    -------
    TYPE
        DESCRIPTION.

    """
    X1, Y1, V1, age1 = streamline(x, y, vx, vy, xs, ys, stepsize=50.0, minvel=mv, timesize=None, maxpts=maxpts_down)
    X2, Y2, V2, age2 = streamline(x, y, -vx, -vy, xs, ys, stepsize=50.0, minvel=mv, timesize=None, maxpts=maxpts_up)
    X = np.hstack((np.flip(np.asarray(X2))[:-1], np.asarray(X1)))
    Y = np.hstack((np.flip(np.asarray(Y2))[:-1], np.asarray(Y1)))
    age = np.hstack((np.flip(np.asarray(age2))[:-1], np.asarray(age1) + age2[-1]))
    rec = np.c_[X[::-1], Y[::-1]]
    dist_center = [0.0]
    for i in range(1, len(X)):
        dist_center.append(dist_center[i - 1] + np.sqrt((rec[i][0] - rec[i - 1][0]) ** 2 + (rec[i][1] - rec[i - 1][1]) ** 2))
    # if not os.path.exists('../../../velocity/centerline_stor.shp'):
    # need to create this tiff to put in plots
    # glib.dict2shp_pts(data_dir + '/negis/velocity/centerline_eastgrip.shp',
    # {'coords': rec, 'dists': np.array(dist_center),
    #  'ages': np.arange(len(dist_center) - 1, -1, -1)})
    return dist_center, rec, np.flip(age, axis=-1), len(X2)


def get_transect(spacing=1500.0):
    """Get the model domains."""

    tsect = glib.shp2dict("domec_mesh/domec_transect_v2.shp")[0]
    cl = Centerline(tsect)
    cl.firstonly()
    cl.distances(target_spacing=spacing / 100.0)
    ind_edc = np.argmin(np.sqrt((cl.coords[:, 0] - x_edc) ** 2.0 + (cl.coords[:, 1] - y_edc) ** 2.0))
    cl.distance = cl.distance - cl.distance[ind_edc]
    cl.coords = cl.coords[np.abs(cl.distance) <= 30100]
    cl.distance = cl.distance[np.abs(cl.distance) <= 30100]
    ind_edc = np.argmin(np.sqrt((cl.coords[:, 0] - x_edc) ** 2.0 + (cl.coords[:, 1] - y_edc) ** 2.0))

    surf_fn = data_dir + "/antarctica_general/REMA_100m_domec.tif"
    surf = xr.open_dataset(surf_fn, engine="rasterio")
    smooth_surf = surf.map(lambda x: gaussian_filter(x, 25))
    dx = smooth_surf.differentiate("x")
    dy = smooth_surf.differentiate("y")

    x = dx.x.values
    y = dx.y.values
    vx = -dx.band_data.values[0, :, :] * 1000.0
    vy = -dy.band_data.values[0, :, :] * 1000.0

    # Now do the calculation using the averages across the tubes
    # We need the tangents at each point on the FLS--use a center difference
    # for most points but we need to do a left difference and right difference
    # at the ends
    line = cl.coords
    slope = np.hstack(((line[1, 1] - line[0, 1]) / (line[1, 0] - line[0, 0]), (line[2::, 1] - line[:-2:, 1]) / (line[2::, 0] - line[:-2:, 0]), (line[-2, 1] - line[-1, 1]) / (line[-2, 0] - line[-1, 0])))
    slope = -1.0 / slope

    # x_offset should be positive if
    ftw = 800.0
    x_sign = -np.hstack(([-1.0], np.sign(line[2:, 1] - line[:-2, 1]), [1.0]))
    x_sign[x_sign == 0.0] = 1.0
    x_off = x_sign * np.sqrt(ftw**2.0 / (slope**2.0 + 1.0))
    y_off = x_sign * slope * np.sqrt(ftw**2.0 / (slope**2.0 + 1.0))

    # Two edge cases
    x_off[slope == 0.0] = ftw
    y_off[slope == 0.0] = 0.0
    x_off[np.isinf(slope)] = 0.0
    y_off[np.isinf(slope)] = ftw
    x_off[np.isnan(slope)] = 0.0
    y_off[np.isnan(slope)] = ftw

    offset = np.vstack((x_off, y_off)).T
    print("Obtained slope and offset")

    center_east = cl.coords[ind_edc - 200, :]
    center_slope = (cl.coords[ind_edc + 1, 1] - cl.coords[ind_edc - 1, 1]) / (cl.coords[ind_edc + 1, 0] - cl.coords[ind_edc - 1, 0])
    center_off_y = 1.0
    center_off_x = -center_slope * center_off_y
    mag = (center_off_x**2.0 + center_off_y**2.0) ** (1.0 / 2.0)
    scale = 2000.0 / mag
    center_off_x = scale * center_off_x
    center_off_y = scale * center_off_y

    ldist, lcenterline, lage, _ = make_centerlines(x, y, vx, vy, cl.coords[ind_edc, 0] + center_off_x, cl.coords[ind_edc, 1] + center_off_y)

    rdist, rcenterline, rage, _ = make_centerlines(x, y, vx, vy, cl.coords[ind_edc, 0] - center_off_x, cl.coords[ind_edc, 1] - center_off_y)

    ldist, lcenterline2, lage, _ = make_centerlines(x, y, vx, vy, center_east[0] + center_off_x, center_east[1] + center_off_y)

    rdist, rcenterline2, rage, _ = make_centerlines(x, y, vx, vy, center_east[0] - center_off_x, center_east[1] - center_off_y)

    lcenterline = np.vstack((lcenterline2, np.flipud(lcenterline)))
    rcenterline = np.vstack((rcenterline2, np.flipud(rcenterline)))

    if debug:
        plt.figure()
        plt.imshow(surf.band_data[0, :, :], extent=[surf.x[0], surf.x[-1], surf.y[-1], surf.y[0]])
        ll = plt.plot(lcenterline[:, 0], lcenterline[:, 1], linestyle="dotted")
        plt.plot(cl.coords[:, 0], cl.coords[:, 1])
        lr = plt.plot(rcenterline[:, 0], rcenterline[:, 1], linestyle="dashed")
        plt.plot(x_edc, y_edc, marker="d")
        plt.plot(center_east[0], center_east[1], marker="*")
        plt.plot(x_edc - center_off_x, y_edc - center_off_y, marker="s")
        plt.plot(x_edc + center_off_x, y_edc + center_off_y, marker="^")

    print("Neighboring streamlines found")
    l_line = shapely.geometry.LineString(lcenterline)
    r_line = shapely.geometry.LineString(rcenterline)
    if debug:
        print("Lines made")

    l_ends = offset * 10.0 + cl.coords
    r_ends = -offset * 10.0 + cl.coords
    if debug:
        plt.plot(l_ends[:, 0], l_ends[:, 1])
        plt.plot(r_ends[:, 0], r_ends[:, 1])
        print("Intersections found")
        plt.axis("equal")

    int_lines = [shapely.geometry.LineString([l_e, r_e]) for l_e, r_e in zip(l_ends, r_ends)]
    lints = [l_line.intersection(int_line) for int_line in int_lines]
    left_point = line - offset
    right_point = line + offset
    left_point = np.array([intl.coords[0] if isinstance(intl, shapely.geometry.Point) else left_point[j, :] for j, intl in enumerate(lints)])

    rints = [r_line.intersection(int_line) for int_line in int_lines]
    right_point = np.array([intr.coords[0] if isinstance(intr, shapely.geometry.Point) else right_point[j, :] for j, intr in enumerate(rints)])
    if debug:
        plt.plot(left_point[:, 0], left_point[:, 1], color=ll[0].get_color())
        plt.plot(right_point[:, 0], right_point[:, 1], color=lr[0].get_color())
        print("Intersections found")
        plt.axis("equal")

    ftw = gaussian_filter1d(np.sqrt(np.sum((left_point - right_point) ** 2.0, axis=1)), 50)

    if debug:
        plt.figure()
        plt.plot(cl.distance, ftw)
        plt.show()

    bed_fn = data_dir + "/antarctica_general/bedmachine_bed.tif"
    acc_fn = data_dir + "/antarctica_general/snw_final.nc"
    ts_fn = data_dir + "/antarctica_general/ts_mean.nc"

    temp_ds = xr.open_dataset(ts_fn)
    mean_temp = temp_ds["ts"].to_numpy()
    temp_ds.close()

    acc_ds = xr.open_dataset(acc_fn)
    x, y, _ = ll2sps(acc_ds["lat"].to_numpy().flatten(), acc_ds["lon"].to_numpy().flatten())
    acc_interper = LinearNDInterpolator(np.vstack((x, y)).T, acc_ds["snw"].to_numpy()[:, :].flatten())
    acc_ds.close()

    temp_interper = LinearNDInterpolator(np.vstack((x, y)).T, mean_temp.flatten())

    bedint = Interpolator(*glib.gtif2mat_fn(bed_fn, ndv=-9999))
    surfint = Interpolator(*glib.gtif2mat_fn(surf_fn, ndv=-9999))

    line = cl.coords

    bed = bedint.vector_call(line[:, 0], line[:, 1])
    surf = surfint.vector_call(line[:, 0], line[:, 1])
    acc = acc_interper(line[:, 0], line[:, 1])
    acc = acc / acc[ind_edc]
    ts = temp_interper(line[:, 0], line[:, 1])
    ts = ts - ts[ind_edc]

    props = {"coords": line, "bed": bed, "surf": surf, "acc": acc, "ts": ts, "dist": cl.distance, "ftw": ftw, "x_sps": line[:, 0], "y_sps": line[:, 1]}
    for key in props.keys():
        if key not in ["coords", "dist"]:
            # 1D outputs
            np.savetxt("../edc_data/v2/{:s}.txt".format(key), np.vstack((props["dist"], props[key])).T, fmt="%8.5e", header=str(len(props["dist"])))

    glib.dict2shp_pts("domain_v2.shp", props)

    # some things work only on certain machines
    top = np.flipud(np.vstack((cl.distance[::100], V_SCALE * surf[::100])).T)
    bottom = np.vstack((cl.distance[::100], V_SCALE * bed[::100])).T
    fn_base = "domec_v2"
    gmshlib.gmsh_outline(fn_base + ".geo", [bottom[:-1, :], [bottom[-1, :]], top[:-1, :], [top[-1, :]]], RES, cuts=None, cuts_lc=None, points=None, points_lc=None, spline_or_line="Line")

    sp.call(["gmsh", "-1", "-2", fn_base + ".geo", "-algo", "front2d"])
    sp.call(["ElmerGrid", "14", "2", fn_base + ".msh", "-autoclean"])
    # sp.call(['ElmerGrid', '2', '5', fn_base])

    folder = fn_base
    f_in = folder + "/mesh.nodes"
    with open(f_in, "r") as fin:
        lines = fin.readlines()

    with open(f_in, "w") as fout:
        with open(folder + "/meshbackup.nodes", "w") as fbackout:
            for line in lines:
                fbackout.write(line)
                fields = line.split(" ")
                fields[3] = str(float(fields[3]) / V_SCALE)
                new_line = " ".join(fields)
                fout.write(new_line)

    # only works with certain paths, skip
    sp.call(["ElmerGrid", "2", "5", fn_base])


if __name__ == "__main__":
    get_transect()
