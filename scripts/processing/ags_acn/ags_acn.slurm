#!/bin/bash

#SBATCH --job-name=ags_acn
#SBATCH --output=%j.out
#SBATCH --error=%j.err
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=40:00:00
#SBATCH --account=gornish

#source ~/.bashrc
module load anaconda
source ~/.bashrc && conda activate
conda activate vs2

source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/ags_acn/config.sh

export FILE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

