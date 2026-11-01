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
    cd "$SCRIPT_DIR/garfield" || exit 1
    
    git config core.trustctime false
    git config core.checkstat minimal
    
    git checkout master
    git fetch origin master
    git merge FETCH_HEAD --ff-only
    
    cd "$SCRIPT_DIR" || exit 1
fi

# --- PREPARE SCRATCH ---
if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    echo "=== Preparing SCRATCHDIR ==="
    mkdir -p "$SCRATCHDIR/tmp"
    export TMPDIR="$SCRATCHDIR/tmp"
    
    echo "Syncing source..."
    # Exclude build and install directories so we don't unnecessarily copy them as source files
    rsync -a --info=progress2 --exclude='build' --exclude='install' "$SCRIPT_DIR/garfield/" "$SCRATCHDIR/garfield/" || exit 1
    
    if [ -d "$SCRIPT_DIR/garfield/build" ]; then
        echo "Restoring previous build directory to SCRATCHDIR..."
        mkdir -p "$SCRATCHDIR/garfield/build"
        rsync -a --info=progress2 "$SCRIPT_DIR/garfield/build/" "$SCRATCHDIR/garfield/build/"
        
        # Purge hardcoded absolute paths from previous PBS jobs
        find "$SCRATCHDIR/garfield/build" -name "CMakeCache.txt" -delete
        find "$SCRATCHDIR/garfield/build" -name "CMakeFiles" -type d -exec rm -rf {} +
        
        if $DO_BUILD; then
            DO_CONFIGURE=true
        fi
    else
        mkdir -p "$SCRATCHDIR/garfield/build"
    fi
    
    source /etc/profile.d/30_meta_modules.sh
    module add cmake/3.23.1-gcc-10.2.1-gxvea6z
    source /storage/projects/utefx17/SourceCode/ROOT/automatic/root_install/bin/thisroot.sh
    source /storage/projects/utefx17/SourceCode/geant4/automatic/geant4-install/bin/geant4.sh
fi

# --- CONFIGURE ---
if $DO_CONFIGURE; then
    echo "=== Configuring ==="
    cd "$SCRATCHDIR/garfield/build" || exit 1

    LOCAL_QT_CMAKE="$SCRIPT_DIR/../../qt/automatic/qt6/6.11.3/gcc_64/lib/cmake"

    # Replaced $GARFIELD_HOME with ".." to ensure CMake reads source files from the fast scratch disk 
    # instead of traversing back to the network drive if GARFIELD_HOME is set in your environment.
    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$SCRIPT_DIR/garfield/install" \
        -DCMAKE_CXX_STANDARD=23 \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        -DCMAKE_PREFIX_PATH="$LOCAL_QT_CMAKE" \
        .. || exit 2
fi

# --- BUILD ---
if $DO_BUILD; then
    echo "=== Building ==="
    cd "$SCRATCHDIR/garfield/build" || exit 1
    cmake --build . -j "${PBS_NPROCS:-1}" || exit 3
    
    echo "Syncing build directory back to network storage..."
    mkdir -p "$SCRIPT_DIR/garfield/build"
    rsync -a --delete --info=progress2 "$SCRATCHDIR/garfield/build/" "$SCRIPT_DIR/garfield/build/" || exit 3
fi

# --- INSTALL ---
if $DO_INSTALL; then
    echo "=== Installing ==="
    cd "$SCRATCHDIR/garfield/build" || exit 1
    
    if [ -d "$SCRIPT_DIR/garfield/install" ]; then
        rm -rf "$SCRIPT_DIR/garfield/install_old"
        mv "$SCRIPT_DIR/garfield/install" "$SCRIPT_DIR/garfield/install_old"
    fi
    
    cmake --install . || exit 4
    
    echo "Garfield++ installation completed successfully."
fi

# --- CLEANUP ---
if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    echo "=== Cleaning SCRATCHDIR ==="
    rm -rf "$SCRATCHDIR/garfield" "$SCRATCHDIR/tmp"
fi

cd "$START_DIR"
echo "Done."
