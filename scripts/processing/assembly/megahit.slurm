#!/usr/bin/bash
# --------------------------------------------------
# Request resources here
# --------------------------------------------------
#SBATCH --job-name=megahit_array
#SBATCH --output=output_oomsamples_%j.txt
#SBATCH --error=errors_oomsamples_%j.txt
#SBATCH --account=laurameredith
#SBATCH --mail-type=ALL
#SBATCH --mail-user=mtoucedasuarez@hpc.arizona.edu
#SBATCH --partition=standard  
#SBATCH --time=100:00:00
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --cpus-per-task=28
#SBATCH --mem-per-cpu=4000

# --------------------------------------------------
# Load modules here
# --------------------------------------------------
module load anaconda
source ~/.bashrc && conda activate
conda activate megahit

# --------------------------------------------------
# Execute commands here
# --------------------------------------------------

# import config file and read list samples
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/assembly/config.sh
export sample=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

# megahit assembly
megahit -t 10 --presets meta-large -1 $CLEAN_READS_DIR/${sample}_clean_1.fastq -2 $CLEAN_READS_DIR/${sample}_clean_2.fastq -o $OUT_PATH/$sample --out-prefix $sample

# remove intermediate files
rm -R $OUT_PATH/$sample/intermediate_contigs

