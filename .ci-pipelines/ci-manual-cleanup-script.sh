#!/bin/sh

########################################################################
### CI tests for github.com/KineticPreProcessor/KPP                  ###
### NOTE: This script runs CI-tests manually (for testing/debugging) ###
########################################################################

# Get the list of CI test folders in ALL_TESTS
source ./ci-testing-list.sh

# Current directory
this_dir=$(pwd -P)

# Run each test
# Check status of each individual operation and exit if any do not complete
for this_test in $ALL_TESTS; do

    cd ../ci-tests/$this_test
    [ $? -ne 0 ] && exit 1

    echo ""
    echo ">>>>>>>> Making clean for $this_test <<<<<<<<"
    echo ""
    make -f Makefile_$this_test distclean

    # Also remove other output files from the tests
    rm -f *.m
    rm -f Makefile.$this_test

    cd ..

done

# Remove any log files used to store CI-test results
cd $this_dir
rm -fv *log*

# Return w/ success
exit 0
