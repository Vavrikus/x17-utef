#!/bin/bash

# Get the absolute path to the directory this script is located in
SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"

# Go to project root
cd "$SCRIPT_DIR" || exit 1

# Remove build directory
rm -rf build

# Recreate build directory
mkdir build
cd build || exit 1

# Ignore files from build in clang-tidy
echo "Checks: '-*'" > .clang-tidy
# sed -i -e '1i\//NOLINTBEGIN' -e '$a\//NOLINTEND' X17_dict.cxx

# Run CMake with Debug build type
cmake -DCMAKE_BUILD_TYPE=Debug ..

# Build with all threads
make -j$(nproc)