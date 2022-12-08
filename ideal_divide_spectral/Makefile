#
# Makefile
# dlilien, 2021-07-22 12:44
#

.PHONY: mesh

all: input_data.so SIA_BC.so

input_data.so: input_data.f90
	elmerf90 -o input_data.so input_data.f90

SIA_BC.so: SIA_BC.f90
	elmerf90 -o SIA_BC.so SIA_BC.f90

mesh: divide_v2/mesh.nodes


divide_v2/mesh.nodes: divide_v2.msh
	ElmerGrid 14 2 divide_v2.msh -autoclean
	./shrink_tall.py

divide_v2.msh: %.msh: %.geo
	gmsh $< -o $@ -1 -2 -algo front2d
