#hipergator
#export FC = ifort
#export FFLAGS = -O3
#export LINAL =
#export MKL = /apps/compilers/intel/2018/1.163/compilers_and_libraries_2018.1.163/linux/mkl/lib/intel64_lin

#local
export FC = gfortran
export FFLAGS = -O3
export LINAL = -L/Users/jamesthorpe/LAPACK/lapack-3.7.0 -llapack -lblas

#both
export bindir = ../bin

objects = gmh.o input.o project.o coupling.o hcbd.o

all: $(objects)
	cp gmh.o $(bindir)/gmh

hcbd.o : hcbd.f90
	$(FC) $(FFLAGS) $(LINAL) -c hcbd.f90

coupling.o : coupling.f90 hcbd.o
	$(FC) $(FFLAGS) $(LINAL) -c coupling.f90 hcbd.o 

project.o : project.f90
	$(FC) $(FFLAGS) -c project.f90 

input.o : input.f90 project.o
	$(FC) $(FFLAGS) -c input.f90 project.o

gmh.o : gmh.f90 input.o project.o coupling.o hcbd.o
	$(FC) $(FFLAGS) $(LINAL) -o gmh.o gmh.f90 input.o project.o coupling.o hcbd.o

clean:
	rm *.o *.mod

