#!/bin/bash

########################################################################
### C-I test cleanup script for github.com/KineticPreProcessor/KPP   ###
###                                                                  ###
### Use this if you have run the C-I tests manually, and wish to     ###
### remove compiled files (*.o, *.mod, *.exe) from C-I test folders. ###
########################################################################

# Load common variables and functions
. ${KPP_HOME}/.ci-pipelines/ci-common-defs.sh

# Current directory
cwd=$(pwd -P)

# Clean up files in each C-I test folder
for this_test in ${GENERAL_TESTS}; do
    clean_ci_test_folder "${this_test}" "${cwd}"
done

# Remove any log files used to store C-I test results
cd $cwd
rm -fv *log*

# Return w/ success
exit 0
