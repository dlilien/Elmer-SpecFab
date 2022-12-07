# Elmer-SpecFab

This repository has all the code for Lilien et al., 2022 *Simulating higher-order fabric structure in a coupled, anisotropic ice-flow model: application to Dome C*.

## Dependencies

To run the model, a main-branch installation of [SpecFab](https://github.com/nicholasmr/specfab) is required. A version of Elmer/Ice with the Specfab interface is needed as well--hopefully this will eventually make it to the main branch, but currently it can be found [here](https://github.com/dlilien/elmerfem/tree/spectral).

Due to several bugs, there are tight requirements on plotting dependencies.
For plotting, Matplotlib <3.6, >3.1 is required.
Cartopy==0.19.0.post1 is known the only version known to work.
