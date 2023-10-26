#!/bin/bash
cd $(dirname $0)

i="meta"

echo "Connecting to ${i}"

ssh ${i} '
    rm /storage/projects/utefx17/martin/data/ion_map/new_sample/*
    rm /storage/projects/utefx17/martin/data/micro_tracks/new_tracks/*
'