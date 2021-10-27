#!/bin/sh

### CI tests for github.com/KineticPreProcessor/KPP ###

# Build and run the runge-kutta (rk) test
cd /kpp/ci-tests/rk
echo "Generating rk test mechanism with KPP..."
../../bin/kpp rk.kpp
[ $? -ne 0 ] && exit 1
echo "Building the rk test..."
make -f Makefile_rk COMPILER=GFORTRAN
[ $? -ne 0 ] && exit 1
echo "Running the rk test..."
./rk.exe
[ $? -ne 0 ] && exit 1
echo "rk test successful!"
echo ""

# Build and run the sdirk test
cd /kpp/ci-tests/sd
echo "Generating sd test mechanism with KPP..."
../../bin/kpp sd.kpp
[ $? -ne 0 ] && exit 1
echo "Building the sd test ..."
make -f Makefile_sd COMPILER=GFORTRAN
[ $? -ne 0 ] && exit 1
echo "Running the sd test..."
./sd.exe
[ $? -ne 0 ] && exit 1
echo "sd test successful!"
echo ""

# Build and run the small_f90 test
cd /kpp/ci-tests/small_f90
echo "Generating the small_f90 test mechanism with KPP..."
../../bin/kpp small_f90.kpp
[ $? -ne 0 ] && exit 1
echo "Building the small_f90 test..."
make -f Makefile_small_f90 COMPILER=GFORTRAN
[ $? -ne 0 ] && exit 1
echo "Running the small_f90 test..."
./small_f90.exe
[ $? -ne 0 ] && exit 1
echo "small_f90 test successful!"
echo ""
