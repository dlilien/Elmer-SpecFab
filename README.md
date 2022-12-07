# Elmer-SpecFab

This repository has all the code for Lilien et al., 2022 *Simulating higher-order fabric structure in a coupled, anisotropic ice-flow model: application to Dome C*.

## Dependencies

To run the model, a main-branch installation of [SpecFab](https://github.com/nicholasmr/specfab) is required. A version of Elmer/Ice with the Specfab interface is needed as well--hopefully this will eventually make it to the main branch, but currently it can be found [here](https://github.com/dlilien/elmerfem/tree/spectral).

Due to several bugs, there are tight requirements on plotting dependencies.
For plotting, Matplotlib <3.6, >3.1 is required.
Cartopy==0.19.0.post1 is known the only version known to work.

There are also data dependencies. To make figure 3, REMA and Bedmachine Antarctica need to be downloaded to the general_data folder. Projection should be EPSG:3413, and each should be a geotif. To make Figure 4, freely available data from Bazin et al., 2013, are needed from Pangaea (these go in the edc_data folder).
