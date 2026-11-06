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

source /etc/profile.d/30_meta_modules.sh
module add cmake/3.23.1-gcc-10.2.1-gxvea6z

LOCAL_QT_CMAKE="/storage/projects/utefx17/SourceCode/qt/automatic/qt6/6.11.3/gcc_64/lib/cmake"

cmake \
 -DCMAKE_BUILD_TYPE=Release \
 -DCMAKE_PREFIX_PATH=$LOCAL_QT_CMAKE \
 ..

# Build all threads
make -j8 |& tee build.log

cd $START_DIR