#!/bin/sh -l
gfortran weno5consv.f90 -o weno5consv
./weno5consv
#gnuplot "plotOutput.p"
