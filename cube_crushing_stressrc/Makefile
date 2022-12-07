.PHONY: square cube

all: square cube

square: square/mesh.nodes

square/mesh.nodes: square.grd
	ElmerGrid 1 2 square.grd

cube: cube/mesh.nodes

cube/mesh.nodes: cube.grd
	ElmerGrid 1 2 cube.grd

# vim:ft=make
#
