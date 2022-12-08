#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.

import subprocess as sp

temps = ["0", "-30", "-5"]

gammas = {"-5": 5.5e-1, "0": 0.725, "-30": 1.17e-1}
lambdas = {"-5": 2.0e-2, "0": 2.0e-2, "-30": 1.68e-2}

threed = ["UC_ZZ", "UE_ZZ"]
twod = ["CC_ZX", "SS_XZ"]
ugrads = twod + threed

folders = {ugrad: "cube" if ugrad in threed else "square" for ugrad in ugrads}
inflow_boundaries = {"UC_ZZ": [5, 6], "UE_ZZ": [1, 2, 3, 4], "CC_ZX": [1, 3], "SS_XZ": [2, 4]}

Nt = 500
eps_0 = 0.1
dt = (1.0 - 0.01 ** (1.0 / Nt)) / eps_0
tsteps_pure = "  Output Intervals(2) = 1000 1\n  Timestep Intervals(2) = 1 {:d}\n  Timestep Sizes(2) = 0.0000001 {:f}\n".format(Nt, dt)
tsteps_simple = "  Output Intervals(2) = 1000 1\n  Timestep Intervals(2) = 1 {:d}\n  Timestep Sizes(2) = 0.0000001 {:f}\n".format(Nt, 3.0 / Nt / eps_0)

tstep_dict = {"UC_ZZ": tsteps_pure, "UE_ZZ": tsteps_pure, "CC_ZX": tsteps_pure, "SS_XZ": tsteps_simple}


vels = {
    "UC_ZZ": '  AIFlow 1 = Variable Coordinate 1\n    Real MATC "(tx - 0.5) * 5.0e-2"\n  AIFlow 2 = Variable Coordinate 2\n    Real MATC "(tx - 0.5) * 5.0e-2"\n  AIFlow 3 = Variable Coordinate 3\n    Real MATC "(tx - 0.5) * -1.0e-1"\n  AIFlow 4 = Real 0.0\n',
    "UE_ZZ": '  AIFlow 1 = Variable Coordinate 1\n    Real MATC "(tx - 0.5) * -5.0e-2"\n  AIFlow 2 = Variable Coordinate 2\n    Real MATC "(tx - 0.5) * -5.0e-2"\n  AIFlow 3 = Variable Coordinate 3\n    Real MATC "(tx - 0.5) * 1.0e-1"\n  AIFlow 4 = Real 0.0\n',
    "CC_ZX": '  AIFlow 1 = Variable Coordinate 1\n    Real MATC "tx * 1.0e-1"\n  AIFlow 2 = Variable Coordinate 2\n    Real MATC "(tx - 0.5) * -1.0e-1"\n  AIFlow 3 = Real 0.0\n',
    "SS_XZ": '  AIFlow 1 = Variable Coordinate 2\n    Real MATC "(tx - 0.5) * 1.0e-1"\n  AIFlow 2 = Real 0.0\n  AIFlow 3 = Real 0.0\n',
}

latrots = {"DDRX": 0.0, "LATROT": 1.0, "FULL": 1.0}


def prepare_template(ugrad, temp, fname):
    template_fn = "template_cube_crush.sif"
    out_fn = "{:s}_{:s}_{:s}.sif".format(ugrad, fname, temp)

    if fname in ["LATROT", "DDRX"]:
        Λ_0 = 0.0
    else:
        Λ_0 = lambdas[temp]

    if fname == "LATROT":
        Γ_0 = 0.0
    else:
        Γ_0 = gammas[temp]

    if ugrad in threed:
        dim = 3
    else:
        dim = 2

    subs = {
        "TEMP": temp,
        "UGRAD": ugrad,
        "FNAME": fname,
        "GAMMA": Γ_0,
        "DIM": str(dim),
        "AIDIM": str(dim + 1),
        "LAMBDA": Λ_0,
        "LATROT": latrots[fname],
        "FOLDER": folders[ugrad],
        "NBDRS": len(inflow_boundaries[ugrad]),
        "BDRS": " ".join([str(i) for i in inflow_boundaries[ugrad]]),
        "VEL": vels[ugrad],
        "TSTEP_LINES": tstep_dict[ugrad],
    }

    with open(template_fn, "r") as fin:
        with open(out_fn, "w") as fout:
            for line in fin:
                for sub in subs:
                    if "{%s}" % sub in line:
                        line = line.format(**subs)
                        break
                fout.write(line)
    return out_fn


if __name__ == "__main__":
    for ugrad in ugrads:
        for temp in ["-5", "-30", "0"]:  # temps:
            if temp in ["-5"]:
                fnames = ["LATROT", "DDRX", "FULL"]  # , 'LATROT']
            else:
                fnames = ["FULL"]
            for fname in fnames:
                fn = prepare_template(ugrad, temp, fname)
                sp.call(["ElmerSolver", fn])
