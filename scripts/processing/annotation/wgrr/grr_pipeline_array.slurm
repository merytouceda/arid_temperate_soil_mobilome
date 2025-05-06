#!/bin/bash

#SBATCH --job-name=wgrr_pipeline
#SBATCH --output=output_%j.txt
#SBATCH --error=errors_%j.txt
#SBATCH --partition=standard
#SBATCH --ntasks=1
#SBATCH --nodes=1
#SBATCH --mem-per-cpu=4000
#SBATCH --cpus-per-task=40
#SBATCH --time=50:00:00
#SBATCH --account=barberan


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


# Load modules
module load anaconda
source ~/.bashrc && conda activate
#conda activate mmseqs2
#module load python/3.9/3.9.10 
conda activate mmseqs2 || { echo "ERROR: Failed to activate conda environment mmseqs2"; exit 1; }
module load python/3.9/3.9.10 || { echo "ERROR: Failed to load Python module"; exit 1; }

# Initialize config.sh and list of files
source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/urban/annotation/grr/mmseqs_method/per_sample/config.sh
check_error "sourcing config.sh"

# Check if important directories exist
check_dir_exists "$BACT_PROTS_DIR"
check_dir_exists "$PLAS_PROTS_DIR"
check_dir_exists "$MMSEQS_DIR"
check_dir_exists "$JOB_DIR"

# Check if IN_LIST is set and exists
if [ -z "$IN_LIST" ]; then
    echo "ERROR: IN_LIST is not defined in config.sh"
    exit 1
fi
check_file_exists "$IN_LIST"

# Get the file for this array task
export FILE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`
if [ -z "$FILE" ]; then
    echo "ERROR: Failed to get filename from line ${SLURM_ARRAY_TASK_ID} in $IN_LIST"
    exit 1
fi

# Check if input files exist
check_file_exists "$BACT_PROTS_DIR/${FILE}.contigs.fa_pro.fa"
check_file_exists "$PLAS_PROTS_DIR/${FILE}.contigs.fa_pro.fa"




# Initialize config.sh and list of files
#source /groups/egornish/mtoucedasuarez/scripts/jobs/mge/urban/annotation/grr/mmseqs_method/per_sample/config.sh
#export FILE=`head -n +${SLURM_ARRAY_TASK_ID} $IN_LIST | tail -n 1`

start_time=$(date +%s)  # Capture start time
echo "Start time: $start_time"

# CLUSTER - Bacteria 
# Run mmseqs clustering of the protein sequences of plasmids and bacteria to remove duplicates
echo "Ok here we go, running clustering of bacteria"
mmseqs easy-cluster $BACT_PROTS_DIR/${FILE}.contigs.fa_pro.fa $BACT_PROTS_DIR/${FILE}_bact_pro_urban tmp --min-seq-id 0.95 --cov-mode 1 -c 0.7 --threads 48
check_error "mmseqs clustering of bacteria"
check_file_exists "$BACT_PROTS_DIR/${FILE}_bact_pro_urban_rep_seq.fasta"

# CLUSTER - Plasmids
echo "Running clustering of plasmids"
mmseqs easy-cluster $PLAS_PROTS_DIR/${FILE}.contigs.fa_pro.fa $PLAS_PROTS_DIR/${FILE}_plas_pro_urban tmp --min-seq-id 0.95 --cov-mode 1 -c 0.7 --threads 48
check_error "mmseqs clustering of plasmids"
check_file_exists "$PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq.fasta"


# FILTER AND PREPARE
# Filter the sequencs from the bacterial protein set that have been identified as plasmids
echo "Filtering the plasmids out of the bacterial protein files"
python3 $JOB_DIR/filter_fasta.py $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq.fasta $BACT_PROTS_DIR/${FILE}_bact_pro_urban_rep_seq.fasta $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq.fasta
check_error "filtering plasmids from bacterial proteins"
check_file_exists "$BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq.fasta"

# Rename the fasta headers so they read PLASMID_.. and CONTIG_..
echo "...changing the header names ..." 
sed 's/^>/\>CONTIG_/' $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq.fasta > $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq_CONTIG.fasta
check_error "filtering plasmids from bacterial proteins"
check_file_exists "$BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq.fasta"

sed 's/^>/\>PLASMID_/' $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq.fasta > $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq_PLASMID.fasta
check_error "renaming plasmid headers"
check_file_exists "$PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq_PLASMID.fasta"


# Concatenate the two files to create the input file for mmseqs
echo " ...and concatenating the files to get ready for alignment..."
cat $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq_CONTIG.fasta $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq_PLASMID.fasta > $MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fasta
check_error "concatenating fasta files"
check_file_exists "$MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fasta"


# Check file size
file_size=$(stat -c%s "$MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fasta")
if [ "$file_size" -eq 0 ]; then
    echo "ERROR: Concatenated fasta file is empty"
    exit 1
fi



# ALL VS ALL ALIGNMENT
# Run mmseqs all against all
echo "Alright! Starting the all vs. all alginment. Good luck!" 
mmseqs easy-search $MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fasta $MMSEQS_DIR/${FILE}_plasmid_bact_pro_genes_for_mmseqsgrr.fasta $MMSEQS_DIR/${FILE}_plas_bact.m8 $MMSEQS_DIR/tmp -s 7.5 --threads 24 --format-output "query,target,qcov,tcov,fident,evalue,bits"  > $MMSEQS_DIR/${FILE}_plas_all.mmseqs.search.log
check_error "mmseqs all-vs-all alignment"
check_file_exists "$MMSEQS_DIR/${FILE}_plas_bact.m8"


# Check if result file has content
if [ ! -s "$MMSEQS_DIR/${FILE}_plas_bact.m8" ]; then
    echo "ERROR: mmseqs produced an empty result file. Check the log at $MMSEQS_DIR/${FILE}_plas_all.mmseqs.search.log"
    exit 1
fi


# FILTER AND WGRR CALCULATION
# Filter out: 
   # Remove same query-target
   # Keep only the best hit (smallest e-value and highest fident) per query
echo "Now I am just filtering the m8 file"
python3 $JOB_DIR/process_m8.py -f -m $MMSEQS_DIR/${FILE}_plas_bact.m8 -o $MMSEQS_DIR/${FILE}_plas_bact_filtered.txt 
check_error "filtering m8 file"
check_file_exists "$MMSEQS_DIR/${FILE}_plas_bact_filtered.txt"

# Finish filtering and calculate wGRR
   # Keep only mixed hits (plasmid-bact or bact-plasmid, but not plasmid-plasmid or bact-bact)
   # Keep only the best bi-directional hits (you are my best hit and I am yours)
   # Calculate wGRR (https://gitlab.pasteur.fr/jgugliel/wgrr)
echo " and calculating the wGRR"
python3 $JOB_DIR/calculate_grr.py -fm -fb -w -m $MMSEQS_DIR/${FILE}_plas_bact_filtered.txt -a $PLAS_PROTS_DIR/${FILE}_plas_pro_urban_rep_seq_PLASMID.fasta -b $BACT_PROTS_DIR/${FILE}_filtered_bact_pro_urban_rep_seq_CONTIG.fasta -o $MMSEQS_DIR/${FILE}_grr.txt
check_error "calculating wGRR"
check_file_exists "$MMSEQS_DIR/${FILE}_grr.txt"


# Record end time and calculate duration
end_time=$(date +%s)
duration=$((end_time - start_time))
echo "Success! Job completed in $duration seconds"
echo "wGRR output is available at: $MMSEQS_DIR/${FILE}_grr.txt"
