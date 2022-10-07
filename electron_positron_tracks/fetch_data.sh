#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

rsync ${i}:/storage/projects/utefx17/martin/electron_positron_tracks/data/* data/ --progress