#!/bin/bash
#PBS -N ionElectronTracks
#PBS -l select=1:ncpus=1:mem=20gb:scratch_local=10gb
#PBS -l walltime=336:00:00
#PBS -m ae

# The 4 lines above are options for scheduling system. Email notification will be sent when the job aborts (a) or ends (e).

# Make sure we are in the directory of this script.
cd "$(dirname "${BASH_SOURCE[0]}")" || exit 1

# Ensure newly created files are group-writable.
umask 0002

# Check required parameters.
[ -n "$PAR1" ] || { echo "Parameter 1 (max_id) is missing." >&2; exit 5; }
[ -n "$PAR2" ] || { echo "Parameter 2 (id) is missing." >&2; exit 5; }

if [ -n "$PAR3" ]; then
    [ -n "$PAR4" ] || { echo "Parameter 4 (angle_bins) is missing." >&2; exit 5; }
    [ -n "$PAR5" ] || { echo "Parameter 5 (energy_bins) is missing." >&2; exit 5; }
fi

# Directories.
DATADIR="/storage/projects/utefx17/martin/data/micro_tracks/new_tracks"
BUILDDIR="/storage/projects/utefx17/martin/build/simulations/micro_tracks"
MAGDIR="$DATADIR/../../elmag"

# Check scratch directory.
test -n "$SCRATCHDIR" || {
    echo "Variable SCRATCHDIR is not set!" >&2
    exit 1
}

# Record job information.
echo "$PBS_JOBID is running on node $(hostname -f) in scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

# Source ROOT, Geant4 and Garfield++ libraries.
source /storage/projects/utefx17/SourceCode/ROOT/install/bin/thisroot.sh
source /storage/projects/utefx17/SourceCode/geant4/geant4-v11.0.3-install/bin/geant4.sh
source /storage/projects/utefx17/SourceCode/garfield/install/share/Garfield/setupGarfield.sh

# Add GSL to LD_LIBRARY_PATH.
export LD_LIBRARY_PATH="/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/lib:${LD_LIBRARY_PATH:-}"

# Prepare scratch directories.
mkdir -p "$SCRATCHDIR/build/simulations/micro_tracks"
mkdir -p "$SCRATCHDIR/data/elmag"

# Copy executable and input files to scratch.
cp "$BUILDDIR/micro_tracks" "$SCRATCHDIR/build/simulations/micro_tracks/" || {
    echo "Error copying executable!" >&2
    exit 2
}

cp "$MAGDIR/VecE2.txt" "$MAGDIR/VecB2.txt" "$SCRATCHDIR/data/elmag/" || {
    echo "Error copying electric/magnetic field maps!" >&2
    exit 2
}

# Move into scratch directory.
cd "$SCRATCHDIR/build/simulations/micro_tracks" || exit 1

# Run microscopic track simulation.
if [ -z "$PAR3" ]; then
    ./micro_tracks > "tracks${PAR2}.out" 2>&1
else
    ./micro_tracks "$PAR1" "$PAR2" "$PAR3" "$PAR4" "$PAR5" \
        > "tracks${PAR2}.out" 2>&1
fi

# Save simulation exit status immediately.
SIM_EXIT_CODE=$?

# Check the exit status of the simulation.
if [ "$SIM_EXIT_CODE" -ne 0 ]; then
    echo "Calculation failed with exit code $SIM_EXIT_CODE!" >&2

    ERROR_DIR="$DATADIR/error${SIM_EXIT_CODE}"
    mkdir -p "$ERROR_DIR"

    sync_with_group utefx17 \
        "tracks${PAR2}.out" \
        "$ERROR_DIR/" \
        --progress || { echo "Log copying failed!" >&2; exit 4; }

    sync_with_group utefx17 \
        tracks_full1.root \
        "$ERROR_DIR/tracks_full${PAR2}.root" \
        --progress || { echo "Full ROOT file copying failed!" >&2; exit 4; }

    sync_with_group utefx17 \
        tracks_small1.root \
        "$ERROR_DIR/tracks_small${PAR2}.root" \
        --progress || { echo "Small ROOT file copying failed!" >&2; exit 4; }

    exit 3
fi

# Simulation was successful. Copy results to DATADIR.
sync_with_group utefx17 \
    "tracks${PAR2}.out" \
    "$DATADIR/" \
    --progress || {
        echo "Failed to copy simulation log!" >&2
        exit 4
    }

sync_with_group utefx17 \
    tracks_full1.root \
    "$DATADIR/tracks_full${PAR2}.root" \
    --progress || {
        echo "Failed to copy full ROOT output!" >&2
        exit 4
    }

sync_with_group utefx17 \
    tracks_small1.root \
    "$DATADIR/tracks_small${PAR2}.root" \
    --progress || {
        echo "Failed to copy small ROOT output!" >&2
        exit 4
    }

# Clean scratch directory.
clean_scratch
