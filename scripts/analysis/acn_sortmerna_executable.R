# Calculate average copy number (ACN)
# code is modified from: https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# Set up command line arguments
args<-commandArgs(TRUE)

## Input files and directories
wd_path=args[1] # working directory path
blast_file_pattern=args[2] # e.g. "*blast*"
sortmerna_bact_arch_file_path=args[2] # smr_v4.3_bac_archaea.txt

setwd(wd_path)

# load required libraries
library(data.table)

# read blast results of sortmerna
blast_filename <- list.files(pattern = blast_file_pattern)

blast <- lapply(blast_filename, function(i){
  i <- paste("", i, sep = "")
  fread(i,header = F,  sep = "\t")
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
