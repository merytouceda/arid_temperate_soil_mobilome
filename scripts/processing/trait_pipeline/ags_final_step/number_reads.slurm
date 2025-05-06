#!/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=number_reads
#SBATCH --output=number_reads.out
#SBATCH --error=number_reads.err 
#SBATCH --account=gornish
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard
#SBATCH --ntasks=10
#SBATCH --nodes=1
#SBATCH --mem=50gb
#SBATCH --time=02:00:00


cd /xdisk/egornish/mtoucedasuarez/base/trimmomatic/combined

for file in ./*clean_1.fastq ; do

echo "$file: $(awk 'NR % 4 == 1' "$file" | wc -l)" >> number_reads_sample.txt;

done
