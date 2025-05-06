#!/bin/bash -l

# load job configuration
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/mapping_contigs/virus/genomad_sequences/config.sh

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
echo "launching bwa_species.slurm"

JOB_ID=`sbatch --job-name bwa_species -a 1-$NUM_JOB bwa_species.slurm`
