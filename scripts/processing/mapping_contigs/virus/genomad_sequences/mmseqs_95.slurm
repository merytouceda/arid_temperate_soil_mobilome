#!/bin/bash

#SBATCH --job-name=mmseqs2_base_votus
#SBATCH -o %x.out
#SBATCH -e %x.err 
#SBATCH --nodes=1
#SBATCH --ntasks=94
#SBATCH --mem=470gb
#SBATCH --time=10:00:00
#SBATCH --partition=standard
#SBATCH --account=barberan


module load anaconda
source ~/.bashrc && conda activate
conda activate mmseqs2


cd /xdisk/egornish/mtoucedasuarez/base/genomad/combined/all_contigs_virus_sequences/modified_headers

# cluster based on 95% identity and 90% overlap, we use --cov-mode 1 to being related to cd-hit
# See https://mmseqs.com/latest/userguide.pdf
mmseqs easy-cluster /xdisk/egornish/mtoucedasuarez/base/genomad/combined/all_contigs_virus_sequences/modified_headers/all_contigs_virus.fna  base_virus_votus tmp --min-seq-id 0.90 --cov-mode 1 -c 0.8 --threads 94
