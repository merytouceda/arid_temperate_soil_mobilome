#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=bwa_index_votus_genomad
#SBATCH --output=%x.out
#SBATCH --error=%x.err  
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=28
#SBATCH --nodes=1
#SBATCH --mem=168gb
#SBATCH --time=50:00:00

# --------------------------------------------------
# Load modules here
# --------------------------------------------------

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

IN_DIR="/xdisk/egornish/mtoucedasuarez/base/genomad/combined/all_contigs_virus_sequences/modified_headers"
OUT_DIR="/xdisk/egornish/mtoucedasuarez/base/mapping_contigs/combined/virus/genomad_sequences"

# generate index
~/software/bwa/bwa index \
$IN_DIR/base_virus_votus_rep_seq.fasta \
-p $OUT_DIR/votus_genomad_index
