#!/bin/bash

#SBATCH --job-name=hotspot
#SBATCH --output=output_%x.txt
#SBATCH --error=errors_%x.txt
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=40:00:00
#SBATCH --account=gornish

#source ~/.bashrc
#module load anaconda
source ~/.bashrc && conda activate
conda activate hotspot

OUT_DIR="/xdisk/egornish/mtoucedasuarez/base/annotation/hotspot"
IN_DIR="/xdisk/egornish/mtoucedasuarez/base/cluster_ptus"
SOFTWARE_DIR="/groups/egornish/software/HOTSPOT/HOTSPOT"


cd $OUT_DIR

python $SOFTWARE_DIR/preprocessing.py --fasta $IN_DIR/base_plasmid_ptus_rep_seq.fasta --database $SOFTWARE_DIR/database --model_path $SOFTWARE_DIR/models --threads 40

python $SOFTWARE_DIR/hotspot.py --model_path $SOFTWARE_DIR/models --device cpu --threads 40

