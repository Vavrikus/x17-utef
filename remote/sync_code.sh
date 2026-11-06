#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

rsync bashrc.sh ${i}":~/.bashrc" --progress
rsync bash_profile.sh ${i}":~/.bash_profile" --progress 
rsync rebuild_root.sh ${i}":/storage/projects/utefx17/SourceCode/ROOT/automatic" --progress
rsync install_qt.sh ${i}":/storage/projects/utefx17/SourceCode/qt/automatic" --progress
rsync rebuild_geant.sh ${i}":/storage/projects/utefx17/SourceCode/geant4/automatic" --progress
rsync rebuild_garfield.sh ${i}":/storage/projects/utefx17/SourceCode/garfield/automatic" --progress

cd ..
rsync include source simulations reconstruction tests LinkDef.h CMakeLists.txt rebuild_meta.sh ${i}":/storage/projects/utefx17/martin" --info=progress2 --recursive --delete