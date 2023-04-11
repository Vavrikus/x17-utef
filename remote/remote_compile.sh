#!/bin/bash
# cd $(dirname $0)

# i="meta"

# echo "Connecting to ${i}"

# ssh ${i} '
#     source /etc/profile.d/20_meta_modules.sh
#     cd /storage/projects/utefx17/martin/electron_positron_tracks
#     cd build
#     rm -r CMakeCache.txt CMakeFiles/
#     cmake .. -DCMAKE_PREFIX_PATH="/storage/projects/utefx17/SourceCode/libgsl-dev/gsl-2.7-install/;/storage/projects/utefx17/SourceCode/ROOT/install/cmake"
#     make

#     cd ..
#     chmod +x ion_single.sh ion_multi.sh
# '

echo "Script is not up to date."