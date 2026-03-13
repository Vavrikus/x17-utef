#!/bin/bash
#PBS -N clustering
#PBS -l select=1:ncpus=1:mem=1gb:scratch_local=1gb
#PBS -l walltime=00:20:00 
#(((PBS -m ae)))

# The 4 lines above are options for scheduling system. Email notification will be sent when the job aborts (a) or ends (e).

# Ensure necessary parameters are provided
#[ -z $PAR1 ] && echo "Parameter 1 (input_file) is missing." && exit 5
#[ -z $PAR2 ] && echo "Parameter 2 (output_file) is missing." && exit 5

DATADIR=/storage/projects/utefx17/X17data/Hexagon
SRCDIR=/storage/projects/utefx17/SourceCode/TPX3_clustering_linux_mac_v1

# Use the PBS_ARRAY_INDEX to extract the specific file for this job
FILE_LIST="$DATADIR/file_list.dat"

# Ensure the array index is set
if [ -z "$PBS_ARRAY_INDEX" ]; then
    echo >&2 "Error: PBS_ARRAY_INDEX is not set. This script must be run as a job array."
    exit 5
fi

# Extract the Nth line from the file list
PAR1=$(sed -n "${PBS_ARRAY_INDEX}p" "$FILE_LIST")

# Generate the output filename
PAR2="${PAR1%.txt}.root"

# Validate that we actually grabbed a file
[ -z "$PAR1" ] && echo >&2 "Error: Could not read file at index $PBS_ARRAY_INDEX" && exit 5

# Append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory.
# This information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually.
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

# Sourcing ROOT.
source /storage/projects/utefx17/SourceCode/ROOT/root_v6.36/root_install/bin/thisroot.sh

# Set the LD_LIBRARY_PATH to include the path to the GSL library.
export LD_LIBRARY_PATH="/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib:$LD_LIBRARY_PATH"

# Test if scratch directory is set. If scratch directory is not set, issue error message and exit.
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# Set a trap to ALWAYS clean the scratch directory, even if the job exits via an error code
trap 'clean_scratch' TERM EXIT

# --- FIX: Move into the scratch directory BEFORE copying files ---
cd $SCRATCHDIR || { echo >&2 "Failed to enter SCRATCHDIR!"; exit 1; }

# Copy input files to the current directory (which is now $SCRATCHDIR)
cp $SRCDIR/TPX3_clusterer_dev . || { echo >&2 "Error while copying executable!"; exit 2; }
cp $DATADIR/raw/$PAR1 . || { echo >&2 "Error while copying input file(s)!"; exit 2; }

./TPX3_clusterer_dev -i $PAR1 -o $PAR2 -t 200 --time_offset 0 --file_format ASCII_KATHRINE

# Check the exit status of the simulation
if [ $? -ne 0 ]; then
    # Simulation ended with an error.
    echo >&2 "Calculation ended up erroneously (with a code $?) !!"
    # Create a subdirectory in DATADIR with a name based on the number of error.
    error_directory="$DATADIR/error$?"
    mkdir -p "$error_directory"

    # Move the output files to the error directory
    cp $PAR2 $error_directory/ || { echo >&2 "Result file copying failed (with a code $?) !!"; }
    exit 3
else
    # Simulation was successful.
    # Move the output files to user's DATADIR or exit in case of failure.
    cp $PAR2 $DATADIR/root/ || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }
fi

# Clean the SCRATCH directory.
# clean_scratch