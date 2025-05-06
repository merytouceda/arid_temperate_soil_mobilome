#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=acn_r
#SBATCH --output=acn_r.out
#SBATCH --error=acn_r.err 
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --time=02:00:00

module load R

ACN_DIR="/xdisk/egornish/mtoucedasuarez/base/trait/acn"
AGS_DIR="/xdisk/egornish/mtoucedasuarez/base/trait/ags"
JOB_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/trait_pipeline/acn_final_steps"

# Usage: 
# Rscript acn_sortmerna_executable.R [working directory] [pattern blast] [sortmerna bac_arch list]
#Rscript $JOB_DIR/acn_executable.R $ACN_DIR "*coverage*" base

# Rscript acn_sortmerna_executable.R [working directory] [pattern coverage] [AGS file] [name project]
Rscript $JOB_DIR/acn_mine.R $ACN_DIR "*coverage*" $AGS_DIR/base_ags_result.txt  base
