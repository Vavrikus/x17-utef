#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

REMOTE="meta"
PROJECT="/storage/projects/utefx17"

# Apply the same ownership and permission policy as sync_with_group.
RSYNC_GROUP_OPTS=(
    --no-owner
    --no-group
    --chown=:utefx17
    --chmod=ug=rwX
)

echo "Connecting to ${REMOTE}"

# User configuration
rsync -a --progress \
    bashrc.sh "${REMOTE}:~/.bashrc"

rsync -a --progress \
    bash_profile.sh "${REMOTE}:~/.bash_profile"

# Build scripts
rsync -a "${RSYNC_GROUP_OPTS[@]}" --progress \
    rebuild_root.sh \
    "${REMOTE}:${PROJECT}/SourceCode/ROOT/automatic/"

rsync -a "${RSYNC_GROUP_OPTS[@]}" --progress \
    install_qt.sh \
    "${REMOTE}:${PROJECT}/SourceCode/qt/automatic/"

rsync -a "${RSYNC_GROUP_OPTS[@]}" --progress \
    rebuild_geant.sh \
    "${REMOTE}:${PROJECT}/SourceCode/geant4/automatic/"

rsync -a "${RSYNC_GROUP_OPTS[@]}" --progress \
    rebuild_garfield.sh \
    "${REMOTE}:${PROJECT}/SourceCode/garfield/automatic/"

# Project source code
cd ..

rsync -a "${RSYNC_GROUP_OPTS[@]}" \
    --info=progress2 \
    --delete \
    include source simulations reconstruction tests \
    LinkDef.h CMakeLists.txt rebuild_meta.sh \
    "${REMOTE}:${PROJECT}/martin/"

echo "Done."
