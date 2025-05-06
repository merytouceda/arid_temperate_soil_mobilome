# Calculate average genome size (ags)
# code is modified from: https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# Set up command line arguments
args<-commandArgs(TRUE)

## Input files and directories
wd_path=args[1] # working directory path
cog_length_file_path=args[2]
cog_count_file_pattern=args[3] # e.g. "*single_cogs_count.tsv"
total_basepairs_file_pattern=args[4] # e.g. "*totalbp.txt"
name=args[5] # name of the sample to add to the output


# load required libraries
library(data.table)

# set working directory
setwd(wd_path)

# Read data
# ------------------------------------------------------------------------------------------------------------
# 1. Read the length of 35 single-copy genes (one file)
cog.length <- read.table(file = cog_length_file_path, header = T,  sep = "\t", row.names = 1)



# 2. Read annotation results of 35 single-copy genes (multiple files)
sample.cog.filename <- list.files(pattern = cog_count_file_pattern)

sample.cog <- lapply(sample.cog.filename, function(i){
  i <- paste("", i, sep = "")
  read.table(i,header = T,  sep = "\t", row.names = 1)
  })



# 3. Read results of total number of base pairs (multiple files)
sample.totalbp.filename <- list.files(pattern = total_basepairs_file_pattern)

sample.totalbp <- lapply(sample.totalbp.filename, function(i){
  i <- paste("", i, sep = "")
  read.table(i,header = F)
  })

# convert list to data.table
sample.totalbp <- rbindlist(sample.totalbp)

# convert data.table to data.frame
sample.totalbp <- as.data.frame(sample.totalbp)

# rename rows and columns
sample.name <- gsub("_totalbp.txt", "", sample.totalbp.filename)

rownames(sample.totalbp) <- sample.name

colnames(sample.totalbp) <- "totalbp"

names(sample.cog) <- sample.name


# Calculations
# ------------------------------------------------------------------------------------------------------------

# 1. Calculate the coverage of the 35 single-copy genes, copy from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# in case if coverage of some single-copy genes is zero
for (i in names(sample.cog)) {
  if (sum(!rownames(cog.length) %in% rownames(sample.cog[[i]])) > 0 ) {
    j <- !rownames(cog.length) %in% rownames(sample.cog[[i]])
    missing <- rownames(cog.length)[j]
    cov <- rep(0,length(missing))
    sample.cog[[i]][missing,] <- cov
  }
}



# 2. Calculate the coverage of the 35 single-copy genes, that is the number of genomes

genome.number <- lapply(sample.cog, function(x) {mean(x[rownames(cog.length), "cov"]/cog.length$value)})

genome.number <- unlist(genome.number)



# 3. Calculate ags
ags.result <- sample.totalbp

ags.result$genome_number <- genome.number

ags.result$ags <- (ags.result$totalbp/ags.result$genome_number)/1000000

# 4. Write the result to file
ags_file_name <- paste(name, "ags_result.txt", sep="_")
write.table(ags.result, file = ags_file_name, quote=F, sep="\t", col.names = NA)
