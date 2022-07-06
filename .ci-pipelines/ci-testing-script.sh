#!/bin/sh

########################################################################
### CI tests for github.com/KineticPreProcessor/KPP                  ###
########################################################################

# Get the list of CI test folders in ALL_TESTS
source /kpp/.ci-pipelines/ci-testing-list.sh

# Run each test
# Check status of each individual operation and exit if any do not complete
for this_test in $ALL_TESTS; do

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

done

# Return w/ success
echo ""
echo ">>>>>>>> All tests finished succesfully! <<<<<<<<"
echo ""
exit 0
