#!/bin/sh

########################################################################
### CI tests for github.com/KineticPreProcessor/KPP                  ###
### NOTE: This script runs CI-tests manually (for testing/debugging) ###
########################################################################

# Get the list of CI test folders in $ALL_TESTS
source ./ci-testing-list.sh

# Check status of each individual operation and exit if any do not complete
for this_test in $ALL_TESTS; do

    cd ../ci-tests/$this_test
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> Generating $this_test mechanism with KPP <<<<<<<<"
    echo ""
    ../../bin/kpp $this_test.kpp
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> Building the $this_test test executable <<<<<<<<<"
    echo ""
    make -j -f Makefile_$this_test COMPILER=GFORTRAN
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> Running the $this_test test <<<<<<<<"
    echo ""
    ./$this_test.exe
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> $this_test test was successful! <<<<<<<<"
    echo ""

    cd ..

done

# Return w/ success
echo ""
echo ">>>>>>>> All tests finished succesfully! <<<<<<<<"
echo ""
exit 0

