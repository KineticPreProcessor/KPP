#!/bin/sh

########################################################################
### CI tests for github.com/KineticPreProcessor/KPP                  ###
########################################################################

# List of tests (add more as necessary; separate each with a space)
all_tests="radau90 rk rktlm ros rosadj rosenbrock90 rostlm saprc2006 sd sdadj small_f90 ros_upcase"

# Run each test
# Check status of each individual operation and exit if any do not complete
for this_test in $all_tests; do

    cd /kpp/ci-tests/$this_test
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> Generating $this_test mechanism with KPP <<<<<<<<"
    echo ""
    ../../bin/kpp $this_test.kpp
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> Building the $this_test test executable <<<<<<<<<"
    echo ""
    make -f Makefile_$this_test COMPILER=GFORTRAN
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> Running the $this_test test <<<<<<<<"
    echo ""
    ./$this_test.exe
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> $this_test test was successful! <<<<<<<<"
    echo ""

done

# Return w/ success
exit 0
