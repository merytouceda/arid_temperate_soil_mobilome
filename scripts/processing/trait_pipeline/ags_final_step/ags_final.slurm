#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=ags_final_vol15
#SBATCH --output=ags_final_vol5.out
#SBATCH --error=ags_final_vol5.err 
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --time=02:00:00

module load R

AGS_DIR="/xdisk/egornish/mtoucedasuarez/base/trait/ags"
JOB_DIR="/groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/trait_pipeline/ags_final_step"
UPROC_DB_DIR="/xdisk/egornish/mtoucedasuarez/databases/uproc"

#Rscript $JOB_DIR/ags_vol15.R -d $AGS_DIR -l $UPROC_DB_DIR/all_cog_lengths.tsv -c "*single_cogs_count.tsv" -b "*totalbp.txt" -n base 
Rscript  $JOB_DIR/ags_vol1.R $AGS_DIR $UPROC_DB_DIR/all_cog_lengths.tsv "*single_cogs_count.tsv" "*totalbp.txt" base
