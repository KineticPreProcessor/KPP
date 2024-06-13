########################################################################
### Common variables and functions for the C-I tests                 ###
########################################################################

#=======================================================================
# Lists of C-I tests
#=======================================================================

# Testing different mechanisms
GENERAL_TESTS="C_rk
C_rosadj
C_sd
C_sdadj
C_small_strato
F90_implicit12
F90_lsode
F90_radau
F90_rk
F90_rktlm
F90_ros
F90_rosadj
F90_ros_autoreduce
F90_rosenbrock
F90_ros_split
F90_rostlm
F90_ros_upcase
F90_saprc_2006
F90_sd
F90_sdadj
F90_seulex
F90_small_strato
"
# Testing if #MINVERSION works
MINVERSION_TEST="X_minver"

# Testing if the Master Chemical Mechanism test works
DO_MCM=1
MCM_TEST="mcm"

#=======================================================================
# Functions
#=======================================================================

# Returns the path to a C-I test folder
function get_ci_test_path() {
    echo "${KPP_HOME}/ci-tests/${1}"
    return
}

# Run a C-I test
function run_ci_test() {

    # Arguments
    this_test=${1}     # Name of test
    return_dir=${2}    # Directory where we will return upon test completion
    extra_cmds=${3}    # Extra commands to pass to compilation

    # Navigate to C-I test folder (or exit if error)
    test_path=$(get_ci_test_path "${this_test}")
    cd $test_path
    [[ $? -ne 0 ]] && exit 1

    echo ""
    echo ">>>>>>>> Generating ${this_test} mechanism with KPP <<<<<<<<"
    echo ""
    ${KPP_HOME}/bin/kpp $this_test.kpp
    [[ $? -ne 0 ]] && exit 1

    echo ""
    echo ">>>>>>>> Building the ${this_test} test executable <<<<<<<<<"
    echo ""
    make -j -f Makefile_$this_test COMPILER=GFORTRAN ${extra_cmds}
    [[ $? -ne 0 ]] && exit 1

    echo ""
    echo ">>>>>>>> Running the ${this_test} test <<<<<<<<"
    echo ""
    ./$this_test.exe
    [[ $? -ne 0 ]] && exit 1

    echo ""
    echo ">>>>>>>> ${this_test} test was successful! <<<<<<<<"
    echo ""

    # Navigate back to original path and return w/ success
    cd $return_dir
}

# Test if the #MINVERSION command works
function run_minversion_ci_test() {

    # Arguments
    this_test=${1}     # Name of test
    return_dir=${2}    # Directory where we will return upon test completion

    # Navigate to test folder (or exit if error)
    test_path=$(get_ci_test_path "${MINVERSION_TEST}")
    cd ${test_path}
    [[ $? -ne 0 ]] && exit 1

    # Build the mechanism with KPP
    echo ""
    echo ">>>>>>>> Generating ${this_test} mechanism with KPP <<<<<<<<"
    echo ""
    ${KPP_HOME}/bin/kpp ${MINVERSION_TEST}.kpp

    # Navigate back to original path and return w/ success
    cd $return_dir
}

# Removes compiled files in C-I test folders
function clean_ci_test_folder()  {

    # Arguments
    this_test=${1}     # Name of test
    return_dir=${2}    # Directory where we will return upon test completion

    # Navigate to test folder (or exit if error)
    test_path=$(get_ci_test_path "${this_test}")
    cd ${test_path}
    [[ $? -ne 0 ]] && exit 1

    echo ""
    echo ">>>>>>>> Making distclean for ${this_test} <<<<<<<<"
    echo ""
    make -f Makefile_${this_test} distclean

    # Also remove other output files from the tests
    rm -f *.m
    rm -f Makefile.${this_test}

    # Navigate back to original path and return w/ success
    cd ${cwd}
}
