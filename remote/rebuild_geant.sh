#!/bin/bash
START_DIR=$(pwd)

# Get the absolute path to the directory this script is located in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

# 1. Parameter Parsing
DO_FETCH=false
DO_CONFIGURE=false
DO_BUILD=false
DO_INSTALL=false

if [ "$#" -eq 0 ]; then
    DO_FETCH=true
    DO_CONFIGURE=true
    DO_BUILD=true
    DO_INSTALL=true
else
    for arg in "$@"; do
        case $arg in
            fetch) DO_FETCH=true ;;
            configure) DO_CONFIGURE=true ;;
            build) DO_BUILD=true ;;
            install) DO_INSTALL=true ;;
            *) echo "Unknown option: $arg"; exit 1 ;;
        esac
    done
fi

if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    if [ -z "$SCRATCHDIR" ]; then
        echo "Error: SCRATCHDIR is not set. The 'configure', 'build', and 'install' steps must be run in a PBS interactive job."
        exit 1
    fi
fi

# --- FETCH ---
if $DO_FETCH; then
    echo "=== Fetching Source ==="
    cd "$SCRIPT_DIR/geant4-source" || exit 1
    
    git config core.trustctime false
    git config core.checkstat minimal
    
    git fetch --tags

    # Find the latest tag that matches standard versioning (e.g., v11.2.1)
    # The regex ensures we ignore -beta or -rc tags. sort -V does a proper semantic version sort.
    LATEST_TAG=$(git tag | grep -E '^v[0-9]+\.[0-9]+\.[0-9]+$' | sort -V | tail -n 1)
    
    if [ -z "$LATEST_TAG" ]; then
        echo "Error: Could not determine latest Geant4 tag."
        exit 1
    fi

    echo "Checking out latest stable Geant4 release: $LATEST_TAG"
    git checkout "$LATEST_TAG"
    
    cd "$SCRIPT_DIR" || exit 1
fi

# --- PREPARE SCRATCH ---
if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    echo "=== Preparing SCRATCHDIR ==="
    mkdir -p "$SCRATCHDIR/tmp"
    export TMPDIR="$SCRATCHDIR/tmp"
    
    echo "Syncing source..."
    # --info=progress2 prints overall transfer progress rather than a file-by-file list
    rsync -a --info=progress2 "$SCRIPT_DIR/geant4-source/" "$SCRATCHDIR/geant4-source/" || exit 1
    
    if [ -d "$SCRIPT_DIR/geant4-build" ]; then
        echo "Restoring previous geant4-build to SCRATCHDIR..."
        mkdir -p "$SCRATCHDIR/geant4-build"
        rsync -a --info=progress2 "$SCRIPT_DIR/geant4-build/" "$SCRATCHDIR/geant4-build/"
        
        # Purge hardcoded absolute paths from previous PBS jobs
        find "$SCRATCHDIR/geant4-build" -name "CMakeCache.txt" -delete
        find "$SCRATCHDIR/geant4-build" -name "CMakeFiles" -type d -exec rm -rf {} +
        
        if $DO_BUILD; then
            DO_CONFIGURE=true
        fi
    else
        mkdir -p "$SCRATCHDIR/geant4-build"
    fi
    
    source /etc/profile.d/30_meta_modules.sh
    # Assuming cmake is the only strictly required module, load it as in the ROOT script. 
    # Add others (like qt, xerces-c) here if needed for Geant4 on MetaCentrum.
    module add cmake/3.23.1-gcc-10.2.1-gxvea6z
fi

# --- CONFIGURE ---
if $DO_CONFIGURE; then
    echo "=== Configuring ==="
    cd "$SCRATCHDIR/geant4-build" || exit 1

    # Automatically points to whatever version is symlinked as 'latest'
    LOCAL_QT_CMAKE="$SCRIPT_DIR/../../qt/automatic/qt6/6.11.3/gcc_64/lib/cmake"

    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$SCRIPT_DIR/geant4-install" \
        -DCMAKE_CXX_STANDARD=23 \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DCMAKE_PREFIX_PATH="$LOCAL_QT_CMAKE" \
        -DGEANT4_INSTALL_DATA=ON \
        -DGEANT4_USE_GDML=ON \
        -DGEANT4_USE_QT=ON \
        -DGEANT4_USE_RAYTRACER_X11=ON \
        -DGEANT4_USE_OPENGL_X11=ON \
        ../geant4-source || exit 2
fi

# --- BUILD ---
if $DO_BUILD; then
    echo "=== Building ==="
    cd "$SCRATCHDIR/geant4-build" || exit 1
    cmake --build . -j "${PBS_NCPUS:-1}" || exit 3
    
    echo "Syncing geant4-build back to network storage..."
    mkdir -p "$SCRIPT_DIR/geant4-build"
    rsync -a --delete --info=progress2 "$SCRATCHDIR/geant4-build/" "$SCRIPT_DIR/geant4-build/" || exit 3
fi

# --- INSTALL ---
if $DO_INSTALL; then
    echo "=== Installing ==="
    cd "$SCRATCHDIR/geant4-build" || exit 1
    
    if [ -d "$SCRIPT_DIR/geant4-install" ]; then
        rm -rf "$SCRIPT_DIR/geant4-install_old"
        mv "$SCRIPT_DIR/geant4-install" "$SCRIPT_DIR/geant4-install_old"
    fi
    
    cmake --install . || exit 4
    
    echo "Geant4 installation completed successfully."
fi

# --- CLEANUP ---
if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    echo "=== Cleaning SCRATCHDIR ==="
    rm -rf "$SCRATCHDIR/geant4-source" "$SCRATCHDIR/geant4-build" "$SCRATCHDIR/tmp"
fi

cd "$START_DIR"
echo "Done."
