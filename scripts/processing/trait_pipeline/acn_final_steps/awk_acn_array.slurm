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


# C) Initialize config.sh and list of files
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/trait_pipeline/acn_final_steps/config.sh
check_error "sourcing config.sh"


# D) Check if important directories exist
check_dir_exists "$ACN_DIR"    


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
check_file_exists "$ACN_DIR/${sample}.bac_archaea"


# G) Set start time
start_time=$(date +%s)  # Capture start time
echo "Start time: $start_time"



# ---------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------  
# Calculate coverage

    # filter blast output from sortmerna, copy from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/acn.sh
    # l: minimum length to be used in 16S rRNA filtering (default 30)
    # i: minimum identity to be used in 16S rRNA filtering (default 85)
    # e: e-value to be used in 16S rRNA filtering (default 1e-15)
    # s: 16S rRNA reference length (default E. cloi: 1542bp)

awk -v l=30 -v i=85 \
    -v e=1e-15 -v s=1542 \
      '{
        if ( $11 <= e && $4 >= l && $3 >= i ) {
        n_nuc = $10 -$9 +1;
        n_nuc = sqrt(n_nuc*n_nuc)
        n_nuc_tot = n_nuc + n_nuc_tot
      }
    } END {
      print n_nuc_tot/s
    }' $ACN_DIR/${sample}.bac_archaea > $ACN_DIR/${sample}_16S_coverage.txt


# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Success! Job completed in $duration seconds"
echo "Output is available at: $ACN_DIR" # change
