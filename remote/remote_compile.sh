#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

ssh ${i} '
    export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH"
    source /etc/profile.d/20_meta_modules.sh
    cd /storage/projects/utefx17/martin/
    cd simulations/ion_map/
    mkdir -p build
    cd build
    rm -r *
    cmake .. -DCMAKE_PREFIX_PATH="/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/;/storage/projects/utefx17/SourceCode/ROOT/install/cmake"
    make
    cd ..
    chmod +x ion_single.sh ion_multi.sh
'