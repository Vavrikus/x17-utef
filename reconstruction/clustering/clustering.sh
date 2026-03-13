#!/bin/bash
# Move to the raw data directory
cd /storage/projects/utefx17/X17data/Hexagon/raw/

echo "ERROR: This script is currently unsafe as it may mess with permissions. Stopping..."
exit -1

# Create a master list of all text files to be processed
ls -1 *.txt > ../file_list.dat

# Count how many files are in the list
NUM_FILES=$(wc -l < ../file_list.dat)

if [ "$NUM_FILES" -eq 0 ]; then
    echo "No .txt files found. Exiting."
    exit 0
fi

echo "Found $NUM_FILES files. Submitting job array..."

cd /storage/projects/utefx17/X17data/Hexagon/
mkdir -p jobs && cd jobs

# Submit the job array using -J 1-N
# (Optional: limit concurrent jobs by adding %X, e.g., -J 1-$NUM_FILES%50 to run max 50 at once)
qsub -J 1-$NUM_FILES /storage/projects/utefx17/martin/reconstruction/clustering/clustering_single.sh