.PHONY: run

SPECFABDIR=/home/dlilien/sw/specfab/src

SRC=$(wildcard *_L*.f90)
EXECS=$(patsubst %.f90,%,$(SRC))
RUNNERS=$(patsubst %.f90,%.o,$(SRC))

all: $(EXECS)

$(EXECS): %: %.f90
	gfortran $< -I$(SPECFABDIR) -L$(SPECFABDIR) -lspecfab -O2 -mcmodel=medium -lm -llapack -lblas -lnetcdff -I/usr/include -L/usr/lib -ffree-line-length-none -m64 -Wall -fPIC -Wno-integer-division -Wno-unused-dummy-argument -o $@

run: $(RUNNERS)

$(RUNNERS): %.o: %
	./$<

clean:
	rm -f DDRX_L?
	rm -f DDRX_L??
	rm -f LATROT_L?
	rm -f LATROT_L??
	rm -f FULL_L?
	rm -f FULL_L??
	rm -f solutions_L?/*.nc
	rm -f solutions_L??/*.nc
