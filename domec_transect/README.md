# Modeling of Dome C

This runs the real-data modeling.

## Meshing and input data

The data must be read onto the model domain--to avoid the need to download full datasets, the already-interpolated data for the necessary parameters are included in this repository. To fully redo this work, one need only run `real_mesh.py [data_dir]` after downloading the needed data to data_dir (which can be found by going through the tifs listed in real_mesh.py)

## Auxiliary functions

Necessary libraries to read in these inputs can be compiled with `make all`.

## Simulations

There are two simulations of Dome C contained in the two sifs (one using ice-core calibrated parameters, one using lab calibrated parameters). Each will take weeks to run, and is called with `ElmerSolver`.
