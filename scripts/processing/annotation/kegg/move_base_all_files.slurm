#!/bin/bash -l                                                                                                                                                               
#SBATCH --job-name=move_all_base_hmm                                                                                                                                           
#SBATCH --account=laurameredith
#SBATCH --output=move_base_all_hmm.out                                                                                                                                               
#SBATCH --error=move_base_all_hmm.err                                                                                                                                               
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu                                        
#SBATCH --partition=standard                                                                                                                                                 
#SBATCH --time=5:00:00                                                                                                                                                      
#SBATCH --ntasks=24                                                                                                                                                           
#SBATCH --nodes=1                                                                                                                                                            
#SBATCH --mem=120gb

for i in {081..100}; do
    mv /xdisk/egornish/mtoucedasuarez/databases/kegg/profiles/*_base_base_all_pro_rep_seq.part-${i}.fasta.txt  /xdisk/egornish/mtoucedasuarez/databases/kegg/profiles/base_all_protein_kegg
done
