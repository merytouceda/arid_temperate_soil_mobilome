#!/bin/bash

#SBATCH --job-name=awk_acn
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=50:00:00
#SBATCH --account=barberan

# --------------------------------------------------------------------------------------------------------------------------------- 
# 0. SET UP

# A) Create error check functions
# Error handling function
check_error() {
    local exit_code=$?
    local command=$1
    if [ $exit_code -ne 0 ]; then
        echo "ERROR: $command failed with exit code $exit_code at line ${BASH_LINENO[0]}"
        exit $exit_code
    fi
}

# Function to check if file exists
check_file_exists() {
    local file=$1
    if [ ! -f "$file" ]; then
        echo "ERROR: Required file does not exist: $file"
        exit 1
    fi
}

# Function to check directory
check_dir_exists() {
    local dir=$1
    if [ ! -d "$dir" ]; then
        echo "ERROR: Required directory does not exist: $dir"
        exit 1
    fi
}



# B) Load modules
module load anaconda
source ~/.bashrc && conda activate
conda activate seqkit || { echo "ERROR: Failed to activate conda environment seqkit"; exit 1; } 

# C) Initialize config.sh and list of files
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/trait_pipeline/gc_content/config.sh
check_error "sourcing config.sh"


# D) Check if important directories exist
check_dir_exists "$GC_DIR"
check_dir_exists "$INTERLEAVED_READS_DIR"    


# E) Check if IN_LIST is set and exists
if [ -z "$IN_LIST" ]; then
    echo "ERROR: IN_LIST is not defined in config.sh"
    exit 1
fi
check_file_exists "$IN_LIST"

# Get the file for this array task
export sample=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`
if [ -z "$sample" ]; then
    echo "ERROR: Failed to get filename from line ${SLURM_ARRAY_TASK_ID} in $IN_LIST"
    exit 1
fi


# F) Check if input files exist
check_file_exists "$INTERLEAVED_READS_DIR/${sample}_interleave.fa"


# G) Set start time
start_time=$(date +%s)  # Capture start time
echo "Start time: $start_time"



# ---------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------  
# Calculate coverage

seqkit fx2tab $INTERLEAVED_READS_DIR/${sample}_interleave.fa -g -i -n -j 94 > $GC_DIR/${sample}_gc.txt

# mean GC
awk '{ sum += $2 } END { if (NR > 0) print sum / NR }' $GC_DIR/${sample}_gc.txt > $GC_DIR/${sample}_gc_mean.txt

# GC standard deviation
awk '{sum+=$2; sumsq+=$2*$2} END {print sqrt(sumsq/NR - (sum/NR)**2)}' $GC_DIR/${sample}_gc.txt > $GC_DIR/${sample}_gc_var.txt

#rm trait/gc/${sample}_gc.txt
#rm trait/gc/${sample}_interleave.fa


# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Success! Job completed in $duration seconds"
echo "Output is available at: $GC_DIR" # change
