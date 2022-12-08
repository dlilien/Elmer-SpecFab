# Elmer-SpecFab

This repository has all the code for Lilien et al., 2022 *Simulating higher-order fabric structure in a coupled, anisotropic ice-flow model: application to Dome C*.

## Dependencies
The model has been tested on Linux (Ubuntu and CentOS) and on MacOS. It may work on Windows, but it has not been tried.

### Software
Before starting installation, you need an up-to-date Fortran compiler. It needs to support large-rank arrays, so assuming that you use GCC, version 9 is the minimum. We have had success compiling with ifort as well.

To run the model, a main-branch installation of [SpecFab](https://github.com/nicholasmr/specfab) is required. The shared (dynamic) library version must be compiled--this is done by default if you just use `make all` for the build.

A version of Elmer/Ice with the Specfab interface is needed as well--hopefully this will eventually make it to the main branch, but currently it can be found [here](https://github.com/dlilien/elmerfem/tree/spectral). When installing Elmer, you must specify the cmake option -DWITH_SpecFab, and the environment variable SpecFabROOT must be set to the folder containing libspecfab (in addition to the standard glaciology option -DWITH_ElmerIce).

To reproduce meshes, [gmsh](https://gmsh.info) is needed.

My [modeltools](https://github.com/dlilien/modeltools) repository is required for plotting.
Due to several bugs in external libraries, there are tight requirements on the more standard plotting dependencies, too.
For plotting, Matplotlib <3.6, >3.1 is required.
Cartopy==0.19.0.post1 is known the only version known to work.

### Data
There are also data dependencies that are free to download but are not re-hosted here. To make figure 3, REMA and Bedmachine Antarctica need to be downloaded to the general_data folder. Projection should be EPSG:3413, and each should be a geotif. To make Figure 4, freely available data from Bazin et al., 2013, are needed from Pangaea (these go in the edc_data folder). Several plots require the fabric from Durand et al., 2009, to be in the edc_data folder. The pRES plots will require data from Ershadi et al., 2021.

## Reproducing the work

### Calibration
To do all simulations from the ground up, begin with ratefactor_temperature_dependence.py. This will give the temperature dependence of migration recrystallization in the form used in the paper. Further calibration can be done by running plotters_spectral/ratefactor_calibration.py, which will produce Figure 2 of the main text and Supplementary Figure 1. If CALIBRATE is set to True, it will also calibrate the raterfactors--otherwise it will just make the plots.

### Small-scale validation
Next, small-scale validation simulations can be run. These are contained in cube_crushing_stressrc and cube_crushing_strainrc. Each of these folders contains a single python file that will call all the simulations. For comparison, 0d simulations are contained in cube_crushing_stressrc/specfab_0d; the [README in that folder](cube_crushing_stressrc/specfab_0d/README.md) contains information on compilation and running those simulations. After those three sets of simulations are run, Figures 8 and 9 from the paper can be made, to ensure that simulations work as intended.

### Large-scale validation
At this point, everything should be set for large-scale modeling. It is simplest to begin with an idealized case, to confirm everything is working before dealing with real data. An idealized divide is found in ideal_divide_spectral. Further instructions are found in the [README in that folder.](ideal_divide_spectral/README.md)

### Real-data simulations
Finally, the flow-tube simulations can be run. Further instructions are [here.](domec_transect/README.md)
