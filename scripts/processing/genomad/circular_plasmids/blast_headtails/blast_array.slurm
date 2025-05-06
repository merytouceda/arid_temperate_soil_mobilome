#!/bin/bash

#SBATCH --job-name=blast_array
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=1:00:00
#SBATCH --account=barberan


# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load blast/2.13.0 


source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/genomad/circular_plasmids/blast_headtails/config.sh

export SAMPLE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

# blast dbs
makeblastdb -in $IN_DIR/${SAMPLE}_heads.fa -out $DB_DIR/${SAMPLE}_heads_db -dbtype nucl


# for loop over samples:
    
blastn -db $DB_DIR/${SAMPLE}_heads_db -query $IN_DIR/${SAMPLE}_tails.fa -outfmt 6 -out $OUT_DIR/${SAMPLE}_blast_heads_tails.out -max_hsps 1 -evalue 0.0000000001
