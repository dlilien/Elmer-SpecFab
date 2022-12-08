Anisotropic divide
==================

Mesh
----
The mesh is a divide with ice thickness H=2020, alpha=1.0e-3 on either side, and length about 10H on either side of the divide. Bedrock is at y=0. The mesh can be created using `make mesh`; note that there should be a .msh file in the repository, so you should not need gmsh here.

BCs
---

I'm going to use the SIA to get the boundary conditions for the sides; we can do this extra easily since y=0. The files needed to impose this are compiled with `make all`

Running
-------
We want to run three versions:
- Isotropic
- Tensorial fabric
- Spectral fabric

There are three sifs for these three cases, and after the mesh is made and the BCs are compiled, each is just run using `ElmerSolver [sif_fn]`. Note that these will take a long time to run (probably weeks).
