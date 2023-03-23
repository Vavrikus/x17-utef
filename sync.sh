#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

ssh ${i} '
    cd /storage/projects/utefx17/martin
    mkdir -p electron_positron_tracks/build
    mkdir -p mag_data
    mkdir -p electron_positron_tracks/data
'

rsync NSplines.h VectorField.h ${i}":/storage/projects/utefx17/martin" --progress
rsync bashrc.sh ${i}":~/.bashrc" --progress
rsync bash_profile.sh ${i}":~/.bash_profile" --progress
cd electron_positron_tracks
rsync CMakeLists.txt ion_electrons.cpp make_track.cpp reco_track.cpp ion_single.sh ion_multi.sh ${i}":/storage/projects/utefx17/martin/electron_positron_tracks" --progress
cd ../../mag_data
rsync VecE.txt VecB.txt VecE2.txt VecB2.txt ${i}":/storage/projects/utefx17/martin/mag_data" --progress