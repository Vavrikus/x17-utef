#!/bin/bash
cd $(dirname $0) #Makes sure you are in the directory of this script.

[ -z $1 ] && echo "Parameter 1 (max_id) is missing." && exit 6

if [ -z $2 ]; then
   echo "Running random track simulation."
   for i in $(seq 1 $1)
   do
      qsub -v PAR1=$1,PAR2=$i track_single.sh
   done
else
   echo "Parameter 2 (iterations) provided. Running grid-like track simulation."
   [ -z $3 ] && echo "Parameter 3 (angle_bins) is missing." && exit 7
   [ -z $4 ] && echo "Parameter 4 (energy_bins) is missing." && exit 7
   for i in $(seq 1 $1)
   do
      qsub -v PAR1=$1,PAR2=$i,PAR3=$2,PAR4=$3,PAR5=$4 track_single.sh
   done
fi