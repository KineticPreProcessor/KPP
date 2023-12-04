#!/bin/bash

#EOC
#------------------------------------------------------------------------------
#                        The Kinetic PreProcessor (KPP)                       !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: changeVersionNumbers.sh
#
# !DESCRIPTION: Bash script to change the version numbers in the appropriate
#  files in the KPP directory structure.  Run this before releasing a new
#  KPP version.
#\\
#\\
# !CALLING SEQUENCE:
#  $ ./changeVersionNumbers.sh X.Y.Z        # X.Y.Z = KPP version number
#EOP
#------------------------------------------------------------------------------
#BOC

function replace() {
    #========================================================================
    # Replacement for `sed -i -e` that works on both MacOS and Linux
    #
    # 1st argument = regular expression
    # 2nd argument = file to be edited
    #========================================================================
    regex="${1}"
    file="${2}"
    if [[ "x$(uname -s)" == "xDarwin" ]]; then
        sed -i '' -e "${regex}" "${file}"          # MacOS/Darwin
    else
        sed -i -e "${regex}" "${file}"             # GNU/Linux
    fi
}


function exitWithError() {

    #========================================================================
    # Display and error message and exit
    #========================================================================

    echo "Could not update version numbers in ${1}... Exiting!"
    exit 1
}


function main() {

    #========================================================================
    # Replaces the version number in the files listed.
    #
    # 1st argument: New version number to use
    #========================================================================

    # New version number
    version="${1}"

    # Save this directory path and change to root directory
    thisDir=$(pwd -P)
    cd ..

    #========================================================================
    # Update version numbers in various files
    #========================================================================

    # Pattern to match: X.Y.Z
    pattern='[0-9][0-9]*\.[0-9][0-9]*\.[0-9][0-9]*'

    # List of files to replace
    files=(                   \
        "src/gdata.h"         \
        "docs/source/conf.py" \
    )

    # Replace version numbers in files
    for file in ${files[@]}; do
        replace "${pattern}" "${version}" "${file}"
        [[ $? -ne 0 ]] && exitWithError "${file}"
        echo "KPP version updated to ${version} in ${file}"
    done

    #========================================================================
    # Update version number and date in CHANGELOG.md
    #========================================================================

    # Pattern to match: "[Unreleased] - TBD"
    pattern='\[.*Unreleased.*\].*'
    date=$(date -Idate)

    # List of files to replace
    files=(             \
	"CHANGELOG.md"  \
    )

    # Replace version numbers in files
    for file in ${files[@]}; do
	replace "${pattern}" "\[${version}\] - ${date}" "${file}"
        [[ $? -ne 0 ]] && exitWithError "${file}"
        echo "KPP version updated to ${version} in ${file}"
    done

    # Return to the starting directory
    cd "${thisDir}"
}

# ---------------------------------------------------------------------------

# Expect 1 argument, or exit with error
if [[ $# -ne 1 ]]; then
    echo "Usage: ./changeVersionNumbers.sh VERSION"
    exit 1
fi

# Replace version numbers
main "${1}"

# Return status
exit $?
