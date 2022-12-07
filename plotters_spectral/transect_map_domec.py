#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.

"""
Plot Figure 3 (the map)
"""
import yaml
import numpy as np
import matplotlib.pyplot as plt
from modeltools.lib import basemaplib, glib
import matplotlib as mpl

mpl.rcParams["contour.negative_linestyle"] = "solid"

EDC = (1359695.36, -894852.43)
FS = 10

edc_kwargs = {"projection": "stere", "lat_0": -90, "lon_0": 0, "lat_ts": -71, "llcrnrlat": -75.175, "urcrnrlat": -74.86, "llcrnrlon": 124.5, "urcrnrlon": 122.25}


class EDCBasemap(basemaplib.AntarcticaBasemap):
    parallels = np.arange(-80, -70, 1.0)
    meridians = np.arange(100, 150, 1)

    def __init__(self, *args, **kwargs):
        basemaplib.AntarcticaBasemap.__init__(self, *args, **dict(edc_kwargs, **kwargs))

    def label_ll(self, plabels=[False, True, True, True], mlabels=[True, False, False, True], xoff=-0.08, yoff=-0.0, latmax=80.0, color="grey", linecolor="grey", mer_side_offset=None, *args, **kwargs):
        xoffset = xoff * (self.urcrnrx - self.llcrnrx)
        yoffset = yoff * (self.urcrnry - self.llcrnry)
        pl, mr = basemaplib.ProjectedBasemap.label_ll(self, plabels=plabels, mlabels=mlabels, xoffset=xoffset, yoffset=yoffset, color=color, mer_side_offset=mer_side_offset, linecolor=linecolor, *args, **kwargs)
        return pl, mr


def main(width=7.0):
    bm1 = EDCBasemap()
    asp1 = (bm1.llcrnrx - bm1.urcrnrx) / (bm1.llcrnry - bm1.urcrnry)

    bm3 = basemaplib.AntarcticaBasemap()
    asp3 = (bm3.llcrnrx - bm3.urcrnrx) / (bm3.llcrnry - bm3.urcrnry)

    xoff = 0.03
    yoff = 0.06
    yboff = 0.0
    asp1c = (asp1 * (1.0 - yoff)) / (1.0 - xoff)
    height = width / asp1c
    fig, ax = plt.subplots(1, 1, gridspec_kw={"hspace": 0.0, "left": xoff, "right": 1.0, "bottom": 0.0, "top": 1.0 - yoff}, figsize=(width, height))

    bm1 = EDCBasemap(ax=ax)
    bm1.drawparallels(np.arange(-90.0, -60, 0.5), labels=[0, 0, 1, 0], fontsize=FS)
    bm1.drawmeridians(np.arange(110.0, 130.0, 1.0), labels=[0, 0, 1, 0], fontsize=FS)
    ax.text(-0.01, 0.5, "124\N{DEGREE SIGN}E", ha="right", va="center", rotation=90, transform=ax.transAxes, fontsize=10)

    transect_fn = "../domec_transect/domain.shp"
    cm = bm1.plot_tif("../general_data/bedmap2_bed.tif", cmap="BrBG_r", vmin=-300, vmax=300)
    axo = fig.add_axes([xoff, yboff, 0.44, 0.160])
    axo.tick_params(axis="both", which="both", bottom=False, top=False, left=False, right=False, labelbottom=False, labelleft=False)
    cbr = bm1.colorbar(cm, orientation="horizontal", height=0.15, extend="both", ll=(0.25, 0.42 + yboff * height))
    cbr.ax.set_xlabel("Bed elevation (m)", fontsize=10)
    cbr.ax.tick_params(axis="both", which="major", labelsize=10)

    with open("pres_sites.yml") as fin:
        site_data = yaml.load(fin)
    x = np.array([site["x"] for name, site in site_data["sites"].items()])
    y = np.array([site["y"] for name, site in site_data["sites"].items()])
    # mask = np.array([name in site_data['use_sites'] for name in site_data['sites']], dtype=bool)
    mask = np.array([(name in site_data["sparse_sites"] and name != "pREs_Epica") for name in site_data["sites"]], dtype=bool)
    bm1.projected_plot(x[mask], y[mask], color="firebrick", linestyle="none", marker="o", markersize=12, markeredgecolor="k", zorder=11)
    mask_epica = [name == "pREs_Epica" for name in site_data["sites"]]
    mask = np.logical_or(mask, mask_epica)
    # bm1.projected_plot(x[~mask], y[~mask], color='C1', linestyle='none', marker='o', markersize=10, markeredgecolor='k', zorder=11)
    # bm1.projected_plot(x[mask_epica], y[mask_epica], color='C2', linestyle='none', marker='o', markersize=10, markeredgecolor='k', zorder=11)
    for sitename, label in zip(site_data["sparse_sites"], "abcdefghijklmnop"):
        xx = site_data["sites"][sitename]["x"]
        yy = site_data["sites"][sitename]["y"]
        bm1.projected_text(xx, yy, label, fontsize=FS, ha="center", va="center", zorder=12)

    shpdict, shptype = glib.shp2dict(transect_fn)
    x, y = shpdict["coords"][:, 0], shpdict["coords"][:, 1]
    bm1.projected_plot(x, y, color="k", zorder=10)
    bm1.projected_text(x[0] - 100.0, y[0], "A", color="k", fontsize=12, ha="right", va="center")
    bm1.projected_text(x[-1] + 100.0, y[-1], "A'", color="k", fontsize=12, ha="left", va="center")

    ctrs = bm1.contour_tif("../general_data/REMA_100m_dem.tif", colors="0.6", linewidths=0.5, res=1000.0, levels=np.arange(0, 5000, 5.0))
    ax.clabel(ctrs, ctrs.levels, inline=True, fmt="%3.0f m", fontsize=FS)

    # Ice core sites
    bm1.projected_plot(*EDC, color="C2", marker="*", linestyle="none", markersize=12, markeredgecolor="k", zorder=100)
    bm1.projected_text(EDC[0] - 700.0, EDC[1] + 200.0, "EDC (d)", fontsize=FS, ha="right", va="bottom")

    bm1.drawmapscale(lon=123.4, lat=-74.80, length=10, barstyle="fancy", lon0=90, lat0=-71, zorder=99999)

    h = 0.4
    axo = fig.add_axes([xoff, 1 - h - yoff, h * asp3 / asp1c, h])
    bmo = basemaplib.AntarcticaBasemap(ax=axo)
    bmo.drawmapboundary(fill_color="C0")

    bmo.plot_tif("../general_data/bedmap2_surface.tif", vmin=0, vmax=3500, cmap="gray", res=10000.0, zorder=1.5)
    bmo.label_ll(plabels=[0, 0, 0, 0], mlabels=[0, 0, 0, 0], linestyle="solid", linewidth=0.5, latmax=89, color="k", linecolor="k")
    bmo.projected_plot(*EDC, color="C2", marker="*", linestyle="none", markersize=5)
    bmo.projected_text(EDC[0] + 10000.0, EDC[1] - 100000.0, "EDC", fontsize=10, ha="center", va="top", zorder=1000)
    fig.savefig("../plots/figure_3.png", dpi=300)


if __name__ == "__main__":
    main()
