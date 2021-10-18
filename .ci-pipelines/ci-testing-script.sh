#!/bin/sh

# Build and run the small_f90 test
cd /kpp/ci-tests/small.f90
../../bin/kpp small_f90.kpp
make -f Makefile_small_f90 COMPILER=GFORTRAN
./small_f90.exe
