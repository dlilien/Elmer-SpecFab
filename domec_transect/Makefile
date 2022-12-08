#
# Makefile
#

all: SIA_BC.so input_data_v2.so

input_data_v2.so: input_data_v2.f90
	elmerf90 -o input_data_v2.so input_data_v2.f90

SIA_BC.so: SIA_BC.f90
	elmerf90 -o SIA_BC.so SIA_BC.f90
