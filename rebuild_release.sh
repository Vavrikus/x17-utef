#!/bin/bash

START_DIR=$(pwd)

# Get the absolute path to the directory this script is located in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Go to project root
cd "$SCRIPT_DIR" || exit 1

# Remove build directory
rm -rf build

# Recreate build directory
mkdir build
cd build || exit 1

# Run CMake with Debug build type
cmake -DCMAKE_BUILD_TYPE=RelWithDebInfo ..

# Build all threads
make -j$(nproc) |& tee build.log

cd $START_DIR