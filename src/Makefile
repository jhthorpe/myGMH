#hipergator
export FC = ifort
export FFLAGS = -O3
export LINAL =
export MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin
export bindir = ../bin

#both
export pwd = $(shell pwd)
export ML= $(pwd)/linal

objects = gmh.o input.o project.o

all: $(objects)
	cp gmh.o $(bindir)/gmh

project.o : project.f90
	$(FC) $(FFLAGS) -c project.f90 

input.o : input.f90 project.o
	$(FC) $(FFLAGS) -c input.f90 project.o

gmh.o : gmh.f90 input.o project.o
	$(FC) $(FFLAGS) $(LINAL) -o gmh.o gmh.f90 input.o project.o

clean:
	rm *.o *.mod

