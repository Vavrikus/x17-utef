```bash
#!/bin/bash
#PBS -N ionElectronTracks
#PBS -l select=1:ncpus=1:mem=8gb:scratch_local=500mb
#PBS -l walltime=48:00:00
#PBS -m ae

# The 4 lines above are options for scheduling system. Email notification will be sent when the job aborts (a) or ends (e).

# Make sure we are in the directory of this script.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

# Ensure newly created files are group-writable.
umask 0002

# Check all necessary parameters.
[ -n "$PAR1" ] || { echo "Parameter 1 (max_id) is missing." >&2; exit 5; }
[ -n "$PAR2" ] || { echo "Parameter 2 (id) is missing." >&2; exit 5; }
[ -n "$PAR3" ] || { echo "Parameter 3 (step (cm)) is missing." >&2; exit 5; }
[ -n "$PAR4" ] || { echo "Parameter 4 (iterations) is missing." >&2; exit 5; }

# Directories
DATADIR="/storage/projects/utefx17/martin/data/ion_map/new_sample"
BUILDDIR="/storage/projects/utefx17/martin/build/simulations/ion_map"
MAGDIR="$DATADIR/../../elmag"

# Check scratch directory.
test -n "$SCRATCHDIR" || {
    echo "Variable SCRATCHDIR is not set!" >&2
    exit 1
}

# Append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory.
# This information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually.
echo "$PBS_JOBID is running on node $(hostname -f) in scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

# Source ROOT, Geant4 and Garfield++ libraries.
source /storage/projects/utefx17/SourceCode/ROOT/install/bin/thisroot.sh
source /storage/projects/utefx17/SourceCode/geant4/geant4-v11.0.3-install/bin/geant4.sh
source /storage/projects/utefx17/SourceCode/garfield/install/share/Garfield/setupGarfield.sh

# Add GSL to LD_LIBRARY_PATH.
export LD_LIBRARY_PATH="/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib:${LD_LIBRARY_PATH:-}"

# Prepare scratch directories.
mkdir -p "$SCRATCHDIR/build/simulations/ion_map"
mkdir -p "$SCRATCHDIR/data/elmag"

# Copy executable and input files to scratch.
cp "$BUILDDIR/ion_electrons" "$SCRATCHDIR/build/simulations/ion_map/" || {
    echo "Error copying executable!" >&2
    exit 2
}

cp "$MAGDIR/VecE2.txt" "$MAGDIR/VecB2.txt" "$SCRATCHDIR/data/elmag/" || {
    echo "Error copying magnetic/electric field maps!" >&2
    exit 2
}

# Run simulation.
cd "$SCRATCHDIR/build/simulations/ion_map" || exit 1

./ion_electrons "$PAR1" "$PAR2" "$PAR3" "$PAR4" > "ion${PAR2}.out"
status=$?

if [ "$status" -ne 0 ]; then
    echo "Calculation failed with exit code $status!" >&2
    exit 3
fi

# Copy output files to project storage, assigning the correct group
# and permissions to the transferred files.
sync_with_group utefx17 \
    "ion${PAR2}.out" \
    "ion${PAR2}.root" \
    "$DATADIR/" \
    --progress || {
        echo "Result file transfer failed!" >&2
        exit 4
    }

# Clean scratch directory.
clean_scratch
