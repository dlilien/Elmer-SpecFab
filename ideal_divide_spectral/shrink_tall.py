#! /usr/bin/env python
# -*- coding: utf-8 -*-
# vim:fenc=utf-8
#
# Copyright © 2022 David Lilien <david.lilien@umanitoba.ca>
#
# Distributed under terms of the GNU GPL3.0 license.


"""
Make the tall mesh shorter
"""

rescale = 4

folders = ["divide_v2"]

for folder in folders:
    f_in = folder + "/mesh.nodes"
    with open(f_in, "r") as fin:
        lines = fin.readlines()

    with open(f_in, "w") as fout:
        with open(folder + "/meshbackup.nodes", "w") as fbackout:
            for line in lines:
                fbackout.write(line)
                fields = line.split(" ")
                fields[3] = str(float(fields[3]) / rescale)
                new_line = " ".join(fields)
                fout.write(new_line)
