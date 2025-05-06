#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=mobileog_diamond
#SBATCH --output=diamond.out
#SBATCH --error=diamond.err
#SBATCH --account=barberan
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=112gb
#SBATCH --time=125:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load diamond/2.0.9


MOBILE_DB_FA="/xdisk/egornish/mtoucedasuarez/databases/mobileOG-db/mobileOG-db_beatrix-1.6.All.faa"
MOBILE_DB_PATH="/xdisk/egornish/mtoucedasuarez/databases/mobileOG-db"
QUERY_FILE="/xdisk/egornish/mtoucedasuarez/base/prodigal/combined/all_genes_nuc.fa"
OUT_PATH="/xdisk/egornish/mtoucedasuarez/base/mogileOG/combined"

# Make db
diamond makedb --in $MOBILE_DB_FA -d $MOBILE_DB_PATH/mobileOG-diamond-db

# for loop over samples:

diamond blastx --db $MOBILE_DB_PATH/mobileOG-diamond-db  --query $QUERY_FILE  --outfmt 6 --max-target-seqs 1 --evalue 0.0000000001 --out $OUT_PATH/base_combined_mobileOG__121124.out --salltitles  --quiet --more-sensitive


