#!/bin/bash
#PBS -N ionElectronTracks
#PBS -l select=1:ncpus=1:mem=8gb:scratch_local=10gb
#PBS -l walltime=168:00:00 
#(((PBS -m ae)))

# The 4 lines above are options for scheduling system. Email notification will be sent when the job aborts (a) or ends (e).

# Makes sure you are in the directory of this script. Checks all necessary parameters are provided.
cd $(dirname $0)

[ -z $PAR1 ] && echo "Parameter 1 (max_id) is missing." && exit 5
[ -z $PAR2 ] && echo "Parameter 2 (id) is missing." && exit 5
if [ ! -z $PAR3 ]; then
   [ -z $PAR4 ] && echo "Parameter 4 (angle_bins) is missing." && exit 5
   [ -z $PAR5 ] && echo "Parameter 5 (energy_bins) is missing." && exit 5
fi

# Define a DATADIR variable: directory where output will be copied to.
# Define a BUILDDIR and MAGDIR variables: where the input files are taken from. BUILDDIR contains the executable that will be run, MAGDIR electromagnetic data.
DATADIR=/storage/projects/utefx17/martin/data/micro_tracks/new_tracks
BUILDDIR=/storage/projects/utefx17/martin/build/simulations/micro_tracks
MAGDIR=$DATADIR/../../elmag

# Append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory.
# This information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually.
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

# Sourcing ROOT, Geant and Garfield libraries.
source /storage/projects/utefx17/SourceCode/ROOT/install/bin/thisroot.sh
source /storage/projects/utefx17/SourceCode/geant4/geant4-v11.0.3-install/bin/geant4.sh
source /storage/projects/utefx17/SourceCode/garfield/install/share/Garfield/setupGarfield.sh

# Set the LD_LIBRARY_PATH to include the path to the GSL library.
export LD_LIBRARY_PATH="/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib:$LD_LIBRARY_PATH"

# Test if scratch directory is set. If scratch directory is not set, issue error message and exit.
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# Copy all needed input files to scratch directory. If the copy operation fails, issue error message and exit.
mkdir -p $SCRATCHDIR/build/simulations/micro_tracks
mkdir -p $SCRATCHDIR/data/elmag
cp $BUILDDIR/micro_tracks $SCRATCHDIR/build/simulations/micro_tracks || { echo >&2 "Error while copying input file(s)!"; exit 2; }
cp $MAGDIR/VecE2.txt $MAGDIR/VecB2.txt  $SCRATCHDIR/data/elmag || { echo >&2 "Error while copying input file(s)!"; exit 2; }

# Move into scratch directory into the folder that contains the executable.
cd $SCRATCHDIR/build/simulations/micro_tracks

# Run the microscopic track simulation with all necessary parameters and stream its output into a text file. If the calculation ends with an error, issue error message an exit.
if [ -z $PAR3 ]; then
    ./micro_tracks >tracks$PAR2.out || { echo >&2 "Calculation ended up erroneously (with a code $?) !!"; exit 3; }
else
    ./micro_tracks $PAR1 $PAR2 $PAR3 $PAR4 $PAR5 >tracks$PAR2.out || { echo >&2 "Calculation ended up erroneously (with a code $?) !!"; exit 3; }
fi

# Move the output files to user's DATADIR or exit in case of failure.
(cp tracks$PAR2.out $DATADIR/ && cp tracks_full1.root $DATADIR/tracks_full$PAR2.root && cp tracks_small1.root $DATADIR/tracks_small$PAR2.root) || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# Clean the SCRATCH directory.
clean_scratch
