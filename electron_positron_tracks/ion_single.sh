#!/bin/bash
#PBS -N ionElectronTracks
#PBS -l select=1:ncpus=1:mem=4gb:scratch_local=500mb
#PBS -l walltime=1:00:00 
#(((PBS -m ae)))
# The 4 lines above are options for scheduling system: job will run 1 hour at maximum, 1 machine with 4 processors + 4gb RAM memory + 10gb scratch memory are requested, email notification will be sent when the job aborts (a) or ends (e)
cd $(dirname $0) #Makes sure you are in the directory of this script.

[ -z $1 ] && echo "Parameter 1 (max_id) is missing." && exit 5
[ -z $2 ] && echo "Parameter 2 (id) is missing." && exit 5
[ -z $3 ] && echo "Parameter 3 (step (cm)) is missing." && exit 5

# define a DATADIR variable: directory where the input files are taken from and where output will be copied to
DATADIR=/storage/projects/utefx17/martin/electron_positron_tracks/data # substitute username and path to to your real username and path

# append a line to a file "jobs_info.txt" containing the ID of the job, the hostname of node it is run on and the path to a scratch directory
# this information helps to find a scratch directory in case the job fails and you need to remove the scratch directory manually 
echo "$PBS_JOBID is running on node `hostname -f` in a scratch directory $SCRATCHDIR" >> $DATADIR/jobs_info.txt

#loads the Gaussian's application modules, version 03
module add build/ion_electrons

# test if scratch directory is set
# if scratch directory is not set, issue error message and exit
test -n "$SCRATCHDIR" || { echo >&2 "Variable SCRATCHDIR is not set!"; exit 1; }

# copy input file "h2o.com" to scratch directory
# if the copy operation fails, issue error message and exit
#cp $DATADIR/h2o.com  $SCRATCHDIR || { echo >&2 "Error while copying input file(s)!"; exit 2; }

# move into scratch directory
cd $SCRATCHDIR 

# run Gaussian 03 with h2o.com as input and save the results into h2o.out file
# if the calculation ends with an error, issue error message an exit
ion_electrons $1 $2 $3 >ion$2.out || { echo >&2 "Calculation ended up erroneously (with a code $?) !!"; exit 3; }

# move the output to user's DATADIR or exit in case of failure
(cp ion$2.out $DATADIR/ && cp ion$2.root $DATADIR/) || { echo >&2 "Result file(s) copying failed (with a code $?) !!"; exit 4; }

# clean the SCRATCH directory
clean_scratch
