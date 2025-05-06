# Calculate average copy number (ACN)
# code is modified from: https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# Set up command line arguments
args<-commandArgs(TRUE)

## Input files and directories
wd_path=args[1] # working directory path
blast_file_pattern=args[2] # e.g. "*blast*"
sortmerna_bact_arch_file_path=args[3] # smr_v4.3_bac_archaea.txt

# debugging of argument passing
# Diagnostic argument printing
cat("Number of arguments:", length(args), "\n")
cat("Arguments:", args, "\n")
cat("Working Directory:", getwd(), "\n")
cat("BLAST File Pattern:", blast_file_pattern, "\n")
cat("SortMeRNA File Path:", sortmerna_bact_arch_file_path, "\n")

# Set working directory
setwd(wd_path)

# load required libraries
library(data.table)


# 1.  Read blast results of sortmerna
#blast_filename <- list.files(pattern = blast_file_pattern)
blast_filename <- list.files(path = wd_path, pattern = blast_file_pattern, full.names = TRUE) # improved


# Checking that files are being read
print("BLAST files found:")
print(blast_filename)

# Opening files
#blast <- lapply(blast_filename, function(i){
#  i <- paste("", i, sep = "")
#  fread(i,header = F,  sep = "\t")
#  })

# Improved error handling option
blast <- lapply(blast_filename, function(file_path) {
  tryCatch({
    # Verify file exists and is readable
    if (!file.exists(file_path)) {
      stop("File does not exist: ", file_path)
    }
    
    # Attempt to read the file
    fread(file_path, header = FALSE, sep = "\t")
  }, error = function(e) {
    # Detailed error reporting
    cat("Error reading file:", file_path, "\n")
    cat("Error message:", conditionMessage(e), "\n")
    
    # Additional file information
    file_info <- file.info(file_path)
    cat("File size:", file_info$size, "bytes\n")
    cat("File permissions:", file_info$mode, "\n")
    
    # Rethrow the error
    stop(e)
  })
})



# rename list name
sample_name <- gsub(".blast", "", blast_filename)

names(blast) <- sample_name

# read ids of bacteria and archaea
bac_archaea <- read.table(file=sortmerna_bact_arch_file_path, header=T, sep="\t")

# only get bacteria and archaea
blast_bac_archaea <- lapply(blast, function(x){
  x[x$V2%in%bac_archaea$id,]
})

# evalue < 1e-5
blast_bac_archaea <- lapply(blast_bac_archaea, function(x){
  x[x$V11<1e-5,]
})

# export blast results of bacteria and archaea, this is used for calculate 16S coverage
sapply(names(blast_bac_archaea),
       function (x) write.table(blast_bac_archaea[[x]], file=paste(x, "bac_archaea", sep = "."), row.names = F, col.names = F, quote = F, sep = "\t"))
