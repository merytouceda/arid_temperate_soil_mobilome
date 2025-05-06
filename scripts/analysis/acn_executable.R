# Calculate average copy number (ACN)
# code is modified from: https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# Set up command line arguments
args<-commandArgs(TRUE)

## Input files and directories
wd_path=args[1] # working directory path
coverage_file_pattern=args[2] # e.g. "*coverage*"
blast_file_pattern=args[2] # e.g. "*bact_archaea*"
name=args[4] #eg. base

setwd(wd_path)

library(data.table)

sample.coverage.files <- list.files(pattern = coverage_file_pattern)

for (file in sample.coverage.files){
  read.table(file,header = F,  sep = "\t")
}


sample.coverage <- lapply(sample.coverage.filename, function(i){
  i <- paste("", i, sep = "")
  read.table(i,header = F,  sep = "\t")
})

sample.name <- gsub("_16S_coverage.txt", "", sample.coverage.filename)

names(sample.coverage) <- sample.name

sample.coverage <- as.data.frame(rbindlist(sample.coverage))

rownames(sample.coverage) <- sample.name

acn_file_name <- paste(name, "_16S_coverage.txt", sep="_")
write.table(sample.coverage, file=acn_file_name, quote=F, sep="\t", col.names=NA)