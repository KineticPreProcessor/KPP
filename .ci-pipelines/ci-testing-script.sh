#!/bin/sh

### CI tests for github.com/KineticPreProcessor/KPP ###

# # Build and run the runge-kutta (rk) test
# cd /kpp/ci-tests/rk
# ../../bin/kpp rk.kpp
# make -f Makefile_rk COMPILER=GFORTRAN
# ./rk.exe

# # Build and run the sdirk test
# cd /kpp/ci-tests/sd
# ../../bin/kpp sd.kpp
# make -f Makefile_sd COMPILER=GFORTRAN
# ./sd.exe

# Build and run the small_f90 test
cd /kpp/ci-tests/small_f90
../../bin/kpp small_f90.kpp
make -f Makefile_small_f90 COMPILER=GFORTRAN
./small_f90.exe
