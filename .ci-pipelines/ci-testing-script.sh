#!/bin/bash

########################################################################
### C-I test script for github.com/KineticPreProcessor/KPP           ###
########################################################################

# Load common variables and functions
. ${KPP_HOME}/.ci-pipelines/ci-common-defs.sh

# Current directory
cwd=$(pwd -P)

# Run C-I tests with various mechanism + integrator combinations
for this_test in ${GENERAL_TESTS}; do
    run_ci_test "${this_test}" "${cwd}"
done

# Run a C-I test to see if the #MINVERSION command works as advertised
run_minversion_ci_test "${MINVERSION_TEST}" "${cwd}"

# Return w/ success
echo ""
echo ">>>>>>>> All tests finished succesfully! <<<<<<<<"
echo ""

# Return with success
exit 0
