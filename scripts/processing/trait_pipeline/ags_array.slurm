#!/bin/bash

#SBATCH --job-name=ags_pipeline
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
#conda activate mmseqs2
#module load python/3.9/3.9.10 
conda activate sortmerna_env || { echo "ERROR: Failed to activate conda environment mmseqs2"; exit 1; } # change
#module load python/3.9/3.9.10 || { echo "ERROR: Failed to load Python module"; exit 1; } # change


# C) Initialize config.sh and list of files
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/base/trait_pipeline/config.sh
check_error "sourcing config.sh"


# D) Check if important directories exist
check_dir_exists "$CLEAN_READS_DIR"
check_dir_exists "$INTERLEAVED_READS_DIR"    
check_dir_exists "$SORTMERNA_OUTPUT_DIR"
check_dir_exists "$ACN_DIR"    
check_dir_exists "$AGS_DIR"    


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
check_file_exists "$CLEAN_READS_DIR/${sample}_clean_1.fastq"
check_file_exists "$CLEAN_READS_DIR/${sample}_clean_2.fastq"


# G) Set start time
start_time=$(date +%s)  # Capture start time
echo "Start time: $start_time"


# --------------------------------------------------------------------------------------------------------------------------------- 
# 1. Interleave files ~ bbduck 

 # convert R1 and R2 to interleave files
$BBMAP_DIR/bbduk.sh \
in1=$CLEAN_READS_DIR/${sample}_clean_1.fastq \
in2=$CLEAN_READS_DIR/${sample}_clean_2.fastq \
out=$INTERLEAVED_READS_DIR/${sample}_interleave.fa 

check_error "interleaving files"
check_file_exists "$INTERLEAVED_READS_DIR/${sample}_interleave.fa"


# ---------------------------------------------------------------------------------------------------------------------------------                                                                                                                                          
# 2. Run SortmeRNA   

sortmerna \
--ref $SORTMERNA_DB_DIR/smr_v4.3_default_db.fasta \
--reads $INTERLEAVED_READS_DIR/${sample}_interleave.fa \
--workdir $SORTMERNA_OUTPUT_DIR/${sample}_tmp \
--fastx \
-m 3072 \
-e 1e-1 \
--blast 1 \
--num_alignments 1 \
-threads 40

check_error "running sortmerna"

# move blast file to acn
mv $SORTMERNA_OUTPUT_DIR/${sample}_tmp/out/aligned.blast $ACN_DIR/${sample}.blast

# remove intermediate files
rm -R $SORTMERNA_OUTPUT_DIR/${sample}_tmp


# ---------------------------------------------------------------------------------------------------------------------------------
# ---------------------------------------------------------------------------------------------------------------------------------  
# AVERAGE GENOME SIZE (AGS)

# 1. Predict ORF on raw reads ~ FGSpp
$FGSPP_DIR/FGSpp \
-s $INTERLEAVED_READS_DIR/${sample}_interleave.fa \
-o $AGS_DIR/${sample}_fgs \
-w 0 \
-r $FGSPP_DIR/train/ \
-t illumina_5 \
-p 40 \
-c 400

check_error "running FGSPP"

# 2. Calculate the total number of base pairs
egrep -v "^>" $INTERLEAVED_READS_DIR/${sample}_interleave.fa | wc | awk '{print $3-$1}' > $AGS_DIR/${sample}_totalbp.txt
check_error "calculate number of base pairs"


# 3. annotate single copy genes using uproc
$UPROC_DIR/bin/uproc-prot \
    --threads 40 \
    --output $AGS_DIR/${sample}_uproc \
    --preds \
    --pthresh 3 \
    $UPROC_DB_DIR/SINGLE_COPY_COGS_DB \
    $UPROC_DB_DIR/model \
    $AGS_DIR/${sample}_fgs.faa

check_error "running Uproc"

# 4. Get single copy gene count
# This will generate ${sample}_single_cogs_count.tsv, we need this and total number of base pairs to calculate the AGS
cut -f1,3,4,5 -d"," $AGS_DIR/${sample}_uproc | \
    awk 'BEGIN {OFS="\t"; FS=","} {
    if (array_score[$1]) {
    if ($4 > array_score[$1]) {
      array_score[$1] = $4
      array_line[$1, $3] = $2
    }
    } else {
    array_score[$1] = $4
    array_line[$1, $3] = $2
    }
    } END {
    printf "%s\t%s\n", "cog","cov"
    for (combined in array_line) {
    split(combined, separate, SUBSEP)
    array_length[separate[2]]= array_length[separate[2]] + array_line[combined]
    }
    for ( c in array_length ) {
    printf "%s\t%s\n", c,array_length[c]
    }
    }' > $AGS_DIR/${sample}_single_cogs_count.tsv

check_error "running single copy gene count"

    # remove intermediate files
    #rm $INTERLEAVED_READS_DIR/${sample}_interleave.fa 


# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Success! Job completed in $duration seconds"
echo "Output is available at: $AGS_DIR" # change
