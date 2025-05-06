# Calculate average copy number (ACN)
# code is modified from: https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# Set up command line arguments
args<-commandArgs(TRUE)

## Input files and directories
wd_path=args[1] # working directory path
gc_mean_pattern=args[2] # e.g. "*coverage*"
gc_var_pattern=args[3] # e.g. base_ags_result
output_name=args[4] #eg. base

# load necessary libraries
library(tidyverse)

# set working directory
setwd(wd_path)



#------------------------------------------------------------------ GC mean
# read file names
gc_mean_files <- list.files(pattern = gc_mean_pattern)

# initialize empty dataframe
gc_mean_all <- data.frame(matrix(NA, nrow = 1, ncol = 2))
colnames(gc_mean_all) <- c("gc_mean","Sample")

# open all files and read them into dataframe
for (file in gc_mean_files){
  name <- gsub(basename(file), pattern = "_gc_mean.txt", replacement = "")
  d <- read.table(file,header = F,  sep = "\t") %>%
    mutate(Sample = name)
  
  colnames(d) <- c("gc_mean","Sample")
  
  gc_mean_all <- rbind(gc_mean_all, d) %>%
    drop_na()
  
}

#------------------------------------------------------------------ GC var
# read file names
gc_var_files <- list.files(pattern = gc_var_pattern)

# initialize empty dataframe
gc_var_all <- data.frame(matrix(NA, nrow = 1, ncol = 2))
colnames(gc_var_all) <- c("gc_var","Sample")

# open all files and read them into dataframe
for (file in gc_var_files){
  name <- gsub(basename(file), pattern = "_gc_var.txt", replacement = "")
  d <- read.table(file,header = F,  sep = "\t") %>%
    mutate(Sample = name)
  
  colnames(d) <- c("gc_var","Sample")
  
  gc_var_all <- rbind(gc_var_all, d) %>%
    drop_na()
  
}

gc_all <- gc_mean_all %>%
  left_join(.,gc_var_all, by = "Sample")

gc_all_name <- paste(output_name, "gc_result.txt", sep="_")
write.table(gc_all , file=gc_all_name, quote=F, sep="\t", col.names=NA)
