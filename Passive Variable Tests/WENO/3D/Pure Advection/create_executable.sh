#!/bin/sh -l
gfortran -c global.f90
gfortran -c weno5adv.f90 
gfortran global.o weno5adv.o driver.f90 -o driver
./driver
#gnuplot "plotOutput.p"
