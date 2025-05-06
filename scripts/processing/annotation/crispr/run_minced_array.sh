#!/bin/bash -l

# load job configuration
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/annotation/crispr/config.sh

#
#makes sure sample file is in the right place
#
if [[ ! -f "$IN_LIST" ]]; then
    echo "$IN_LIST does not exist. Please provide the path for a list of datasets to process. Job terminated."
    exit 1
fi

# get number of samples to process
export NUM_JOB=$(wc -l < "$IN_LIST")

# submit co_assemblies
echo "launching minced_array.slurm as a job."

JOB_ID=`sbatch --job-name minced_array -a 1-$NUM_JOB minced_array.slurm`
