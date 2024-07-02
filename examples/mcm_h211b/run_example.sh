#!/bin/bash -e

# ============================================================================
# Script to run the MCM minimal example with rosenbrock_h211b_qssa integrator
# ============================================================================

# Export the KPP_HOME folder
export KPP_HOME="../../"

# Remove any leftover files
make -f Makefile_mcm distclean

# Build the MCM mechanism with the rosenbrock_h211b_qssa
$KPP_HOME/bin/kpp mcm.kpp

# Build the executable
gmake -j -f Makefile_mcm EXTERNAL_RATES_F90=constants_mcm.f90

# Run the executable
./mcm.exe

# Display the mixing ratio results
cat mixrat.dat
