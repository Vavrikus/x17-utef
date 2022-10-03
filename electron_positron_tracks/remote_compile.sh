#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

ssh ${i} '
    source /etc/profile.d/20_meta_modules.sh
    cd /storage/projects/utefx17/martin/electron_positron_tracks
    cd build
    cmake ..
    make
'