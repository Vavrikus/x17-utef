#!/bin/bash
START_DIR=$(pwd)

# Get the absolute path to the directory this script is located in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR" || exit 1

PROJECT_GROUP="utefx17"

# 1. Parameter Parsing
DO_FETCH=false
DO_CONFIGURE=false
DO_BUILD=false
DO_INSTALL=false

if [ "$#" -eq 0 ]; then
    # If no arguments are provided, run everything
    DO_FETCH=true
    DO_CONFIGURE=true
    DO_BUILD=true
    DO_INSTALL=true
else
    # Parse specific steps
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

# Ensure heavy operations are only run on compute nodes
if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    if [ -z "$SCRATCHDIR" ]; then
        echo "Error: SCRATCHDIR is not set. The 'configure', 'build', and 'install' steps must be run in a PBS interactive job."
        exit 1
    fi
fi

# --- FETCH ---
# Can be run on the frontend. Network storage operations.
if $DO_FETCH; then
    echo "=== Fetching Source ==="
    cd "$SCRIPT_DIR/root_src" || exit 1

    # NFS Git optimizations: ignore inode changes and minimize stat checks
    git config core.trustctime false
    git config core.checkstat minimal

    git checkout latest-stable || exit 1

    # Fetch and merge instead of pull to avoid some NFS lock overhead
    git fetch origin latest-stable || exit 1
    git merge FETCH_HEAD --ff-only || exit 1

    cd "$SCRIPT_DIR" || exit 1
fi

# --- PREPARE SCRATCH ---
if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    echo "=== Preparing SCRATCHDIR ==="

    # Redirect TMPDIR to avoid 1GB /tmp quota exhaustion during compilation
    mkdir -p "$SCRATCHDIR/tmp"
    export TMPDIR="$SCRATCHDIR/tmp"

    # Sync source to scratch
    rsync -a "$SCRIPT_DIR/root_src/" "$SCRATCHDIR/root_src/" --info=progress2 || exit 1

    # Restore previous build directory for incremental compilation
    if [ -d "$SCRIPT_DIR/root_build" ]; then
        echo "Restoring previous root_build to SCRATCHDIR..."
        mkdir -p "$SCRATCHDIR/root_build"
        rsync -a "$SCRIPT_DIR/root_build/" "$SCRATCHDIR/root_build/" --info=progress2 || exit 1

        # CMake build trees are not relocatable.
        # Purge hardcoded absolute paths from previous PBS jobs.
        find "$SCRATCHDIR/root_build" -name "CMakeCache.txt" -delete
        find "$SCRATCHDIR/root_build" -name "CMakeFiles" -type d -exec rm -rf {} +

        # Reconfiguration is necessary after removing CMake-generated files.
        if $DO_BUILD || $DO_INSTALL; then
            DO_CONFIGURE=true
        fi
    else
        mkdir -p "$SCRATCHDIR/root_build"
    fi

    # Load modules
    source /etc/profile.d/30_meta_modules.sh
    module add cmake/3.23.1-gcc-10.2.1-gxvea6z
    module add libxpm/3.5.12-gcc-10.2.1-czav5y3
    module add gl2ps/1.4.2-gcc-10.2.1-bmc3fym
fi

# --- CONFIGURE ---
if $DO_CONFIGURE; then
    echo "=== Configuring ==="
    cd "$SCRATCHDIR/root_build" || exit 1

    cmake \
        -DCMAKE_BUILD_TYPE=Release \
        -DCMAKE_INSTALL_PREFIX="$SCRIPT_DIR/root_install" \
        -DCMAKE_CXX_STANDARD=23 \
        -DCMAKE_EXPORT_COMPILE_COMMANDS=ON \
        ../root_src || exit 2
fi

# --- BUILD ---
if $DO_BUILD; then
    echo "=== Building ==="
    cd "$SCRATCHDIR/root_build" || exit 1

    cmake --build . -j "${PBS_NCPUS:-1}" || exit 3

    # Sync the build directory back to network storage immediately after a successful build
    echo "Syncing root_build back to network storage for future incremental builds..."
    mkdir -p "$SCRIPT_DIR/root_build"

    # Scratch -> project: ensure correct group ownership and permissions.
    # --delete ensures removed object files don't bloat network storage.
    sync_with_group "$PROJECT_GROUP" \
        "$SCRATCHDIR/root_build/" \
        "$SCRIPT_DIR/root_build/" \
        --delete --info=progress2 || exit 3
fi

# --- INSTALL ---
if $DO_INSTALL; then
    echo "=== Installing ==="
    cd "$SCRATCHDIR/root_build" || exit 1

    # Install locally on scratch first.
    mkdir -p "$SCRATCHDIR/root_install"
    cmake --install . --prefix "$SCRATCHDIR/root_install" || exit 4

    echo "Backing up previous installation..."
    if [ -d "$SCRIPT_DIR/root_install" ]; then
        rm -rf "$SCRIPT_DIR/root_install_old"
        mv "$SCRIPT_DIR/root_install" "$SCRIPT_DIR/root_install_old"
    fi

    echo "Syncing ROOT installation to project directory..."

    # Scratch -> project: ensure correct group ownership and permissions.
    sync_with_group "$PROJECT_GROUP" \
        "$SCRATCHDIR/root_install/" \
        "$SCRIPT_DIR/root_install/" \
        --info=progress2 || exit 4

    echo "ROOT installation completed successfully."
fi

# --- CLEANUP ---
if $DO_CONFIGURE || $DO_BUILD || $DO_INSTALL; then
    echo "=== Cleaning SCRATCHDIR ==="
    rm -rf "$SCRATCHDIR/root_src" \
           "$SCRATCHDIR/root_build" \
           "$SCRATCHDIR/root_install" \
           "$SCRATCHDIR/tmp"
fi

cd "$START_DIR"
echo "Done."
