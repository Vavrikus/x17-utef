#!/bin/bash
cd $(dirname $0) #Makes sure you are in the directory of this script.

[ -z $1 ] && echo "Parameter 1 (max_id) is missing." && exit 6
[ -z $2 ] && echo "Parameter 2 (step (cm)) is missing." && exit 6

for i in $(seq 1 $1)
do
   qsub -- ion_single.sh $1 $i $2
done