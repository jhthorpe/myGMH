#hipergator
export FC = ifort
export FFLAGS = -O3
export LINAL =
export MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
export bindir = /home/james.thorpe/bin/myGMH/bin

#both
export pwd = $(shell pwd)
export ML= $(pwd)/linal

objects = gmh.o input.o

all: $(objects)
	cp gmh.o $(bindir)/gmh

input.o : input.f90
	$(FC) $(FFLAGS) -c input.f90

gmh.o : gmh.f90 input.o
	$(FC) $(FFLAGS) $(LINAL) -o gmh.o gmh.f90 input.o

