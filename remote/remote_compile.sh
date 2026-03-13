#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

ssh ${i} '
    export LD_LIBRARY_PATH="/usr/lib/x86_64-linux-gnu:$LD_LIBRARY_PATH"
    PS1="\[\033[01;32m\]\u@\h\[\033[00m\]:\[\033[01;34m\]\w\[\033[00m\]\$ "
    source /etc/profile.d/30_meta_modules.sh
    module add cmake/3.23.1-gcc-10.2.1-gxvea6z
    source /storage/projects/utefx17/SourceCode/garfield_2603/install/share/Garfield/setupGarfield.sh
    cd /storage/projects/utefx17/martin/
    mkdir -p build
    cd build
    rm -r *
    cmake .. -DCMAKE_PREFIX_PATH="/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/;/storage/projects/utefx17/SourceCode/ROOT/root_v6.36/root_install/cmake"
    make
    cd ../simulations/ion_map/
    chmod +x ion_single.sh ion_multi.sh
    cd ../micro_tracks/
    chmod +x track_single.sh track_multi.sh
    cd ../../reconstruction/clustering/
    chmod +x clustering_single.sh clustering.sh
'