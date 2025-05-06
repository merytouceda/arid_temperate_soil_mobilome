#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=gc_r
#SBATCH --output=gc_r.out
#SBATCH --error=gc_r.err 
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --time=02:00:00

module load R

GC_DIR="/xdisk/egornish/mtoucedasuarez/base/trait/gc"
JOB_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/trait_pipeline/gc_content"

# Usage:
# Rscript gc_combine_samples.R [working directory] [GC mean pattern] [GC var pattern] [name project]
Rscript $JOB_DIR/gc_combine_samples.R $GC_DIR "*_gc_mean.txt" "*_gc_var.txt" base
