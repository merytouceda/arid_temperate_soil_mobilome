#!/bin/bash

#SBATCH --job-name=abricate_mges
#SBATCH --output=%j_mge.out
#SBATCH --error=%j_mge.err
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

source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/annotation/abricate/config.sh

export FILE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

abricate $IN_DIR/${FILE}.contigs.fa_fixednamemobile_contigs.fa > $OUT_DIR/${FILE}_abricate.tab 
