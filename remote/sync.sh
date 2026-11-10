#!/bin/bash
set -e

SCRIPT_DIR="$(cd "$(dirname "${BASH_SOURCE[0]}")" && pwd)"
cd "$SCRIPT_DIR"

./sync_code.sh

i="meta"
REMOTE_DIR="/storage/projects/utefx17/martin/data"

cd ..
rsync -a --progress data/ "${i}:${REMOTE_DIR}/"

ssh "$i" "
    chgrp -R utefx17 '$REMOTE_DIR' &&
    chmod -R u+rwX,g+rwX '$REMOTE_DIR' &&
    find '$REMOTE_DIR' -type d -exec chmod g+s {} +
"
