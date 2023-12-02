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
    run_ci_test "${this_test}" "${cw./cd}" ""
done

# Run the MCM test separately
# NOTE: The MCM test cannot be run on Azure due to memory limitations,
# so test the DO_MCM env var to see if we should run it.
if [[ "$x{DO_MCM}" == "x1" ]]; then
  cmd="EXTERNAL_RATES_F90=constants_mcm.f90 EXTERNAL_RATES_OBJ=constants_mcm.o"
  run_ci_test "${MCM_TEST}" "${cwd}" "${cmd}"
fi

# Run a C-I test to see if the #MINVERSION command works as advertised
run_minversion_ci_test "${MINVERSION_TEST}" "${cwd}"

# Return w/ success
echo ""
echo ">>>>>>>> All tests finished succesfully! <<<<<<<<"
echo ""

# Return with success
exit 0
