#!/bin/bash

#SBATCH --job-name=genomad_base
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=40:00:00
#SBATCH --account=barberan

#source ~/.bashrc
module load anaconda
source ~/.bashrc && conda activate
conda activate mamba
#mamba activate
mamba activate genomad

source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/genomad/config.sh

export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

genomad end-to-end --cleanup $IN_DIR/${SAMPLE}.contigs.fa $OUT_DIR $DATABASE
