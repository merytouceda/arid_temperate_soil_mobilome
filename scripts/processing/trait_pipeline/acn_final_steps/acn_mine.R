
# Calculate average copy number (ACN)
# code is modified from: https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# Set up command line arguments
args<-commandArgs(TRUE)

## Input files and directories
wd_path=args[1] # working directory path
coverage_file_pattern=args[2] # e.g. "*coverage*"
ags_result_file=args[3] # e.g. base_ags_result
output_name=args[4] #eg. base

# load necessary libraries
library(tidyverse)

# set working directory
setwd(wd_path)


#------------------------------------------------------------------Coverage
# read file names
sample_coverage_files <- list.files(pattern = coverage_file_pattern)

# initialize empty dataframe
coverage_all <- data.frame(matrix(NA, nrow = 1, ncol = 2))
colnames(coverage_all) <- c("coverage","Sample")
  
# open all files and read them into dataframe
for (file in sample_coverage_files){
  name <- gsub(basename(file), pattern = "_16S_coverage.txt", replacement = "")
  d <- read.table(file,header = F,  sep = "\t") %>%
    mutate(Sample = name)
    
  colnames(d) <- c("coverage","Sample")
  
  coverage_all <- rbind(coverage_all, d) %>%
    drop_na()

}

# write file
#coverage_file_name <- paste(name, "_16S_coverage.txt", sep="_")
#write.table(coverage_all, file=coverage_file_name, quote=F, sep="\t", col.names=NA)


#------------------------------------------------------------------ACN

# 1. Open ags_result
ags <- read.table(ags_result_file, sep = "\t", header = T)
colnames(ags) <- c("Sample", "totalbp", "genome_number", "ags")

ags <- ags %>%
 mutate(Sample = as.character(as.numeric(Sample)))

acn <- coverage_all %>%
  mutate(Sample = as.character(as.numeric(Sample))) %>%
  left_join(., ags, by = "Sample") %>%
  mutate(acn = coverage/genome_number)


acn_file_name <- paste(output_name, "acn_result.txt", sep="_")
write.table(acn , file=acn_file_name, quote=F, sep="\t", col.names=NA)


