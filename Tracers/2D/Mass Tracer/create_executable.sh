#!/bin/sh
gfortran -c constants.f90
gfortran -fopenmp -c RoeSolver_2D.f90
gfortran -c velocityTracerModule.f90
gfortran -fopenmp constants.o RoeSolver_2D.o velocityTracerModule.o MHD_MH.f90 -o MHD_MH -ffpe-trap=invalid,zero,overflow
