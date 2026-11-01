#!/bin/bash
cd $(dirname $0)
./sync_code.sh

i="meta"

cd ..
rsync data ${i}":/storage/projects/utefx17/martin" --progress --recursive