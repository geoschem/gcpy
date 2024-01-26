#!/bin/bash

#EOC
#------------------------------------------------------------------------------
#                  GEOS-Chem Global Chemical Transport Model                  !
#------------------------------------------------------------------------------
#BOP
#
# !MODULE: changeVersionNumbers.sh
#
# !DESCRIPTION: Bash script to change the version numbers in the appropriate
#  files in the GCPy directory structure.  Run this before releasing a new
#  GCPy version.
#\\
#\\
# !CALLING SEQUENCE:
#  $ ./changeVersionNumbers.sh X.Y.Z        # X.Y.Z = GCClassic version number
#EOP
#------------------------------------------------------------------------------
#BOC

function replace() {

    #========================================================================
    # Function to replace text in a file via sed.
    # 
    # 1st argument: Search pattern
    # 2nd argument: Replacement text
    # 3rd argument: File in which to search and replace
    #========================================================================

    sed -i -e "s/${1}/${2}/" "${3}"
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
    files=(                                                    \
        "docs/source/conf.py"                                  \
	"gcpy/_version.py"                                     \
	"gcpy/benchmark/run_benchmark.py"                      \
	"gcpy/benchmark/modules/run_1yr_fullchem_benchmark.py" \
	"gcpy/benchmark/modules/run_1yr_tt_benchmark.py"       \
    )

    # Replace version numbers in files
    for file in ${files[@]}; do
        replace "${pattern}" "${version}" "${file}"
        [[ $? -ne 0 ]] && exitWithError "${file}"
        echo "GCPy version updated to ${version} in ${file}"
    done

    #========================================================================
    # Update version number and date in CHANGELOG.md
    #========================================================================

    # Pattern to match: "[Unreleased] - TBD"
    pattern='\[.*Unreleased.*\].*'
    date=$(date -Idate)

    # List of files to replace
    file="CHANGELOG.md"

    # Replace version numbers in files
    replace "${pattern}" "\[${version}\] - ${date}" "${file}"
    [[ $? -ne 0 ]] && exitWithError "${file}"
    echo "GCClassic version updated to ${version} in ${file}"

    #========================================================================
    # Update version number in setup.py
    #========================================================================

    # Split version number into an array
    IFS='.' read -r -a version_numbers <<< "${version}"

    # File to replace
    file="setup.py"

    # Replace version numbers
    replace "MAJOR =.." "MAJOR = ${version_numbers[0]}" "${file}"
    [[ $? -ne 0 ]] && exitWithError "${file}"
    replace "MINOR =.." "MINOR = ${version_numbers[1]}" "${file}"
    [[ $? -ne 0 ]] && exitWithError "${file}"
    replace "MICRO =.." "MICRO = ${version_numbers[2]}" "${file}"
    [[ $? -ne 0 ]] && exitWithError "${file}"

    echo "GCPy version updated to ${version} in ${file}"
    
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
