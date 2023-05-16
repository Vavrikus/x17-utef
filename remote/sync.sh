#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

# ssh ${i} '
#     cd /storage/projects/utefx17/martin
#     mkdir -p electron_positron_tracks/build
#     mkdir -p mag_data
#     mkdir -p electron_positron_tracks/data
# '

rsync bashrc.sh ${i}":~/.bashrc" --progress
rsync bash_profile.sh ${i}":~/.bash_profile" --progress 

cd ..
rsync include source simulations data ${i}":/storage/projects/utefx17/martin" --progress --recursive