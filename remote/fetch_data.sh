#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

rsync ${i}:/storage/projects/utefx17/martin/data/ion_map/new_sample/* ../data/ion_map/new_sample/ --progress
rsync ${i}:/storage/projects/utefx17/martin/data/micro_tracks/new_tracks/* ../data/micro_tracks/new_tracks/ --progress