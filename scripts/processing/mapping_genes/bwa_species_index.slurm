#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=bwa_index_all_genes
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

#IN_DIR="/xdisk/egornish/mtoucedasuarez/base/combined_assembly/contigs/plasmid_contigs/circular_10kb"
IN_DIR="/xdisk/egornish/mtoucedasuarez/base/prodigal/combined" # clustered PTUs
OUT_DIR="/xdisk/egornish/mtoucedasuarez/base/mapping_genes/bacteria"

# generate index
~/software/bwa/bwa index \
$IN_DIR/base_all_nuc_rep_seq.fasta \
-p $OUT_DIR/bact_genes_index
