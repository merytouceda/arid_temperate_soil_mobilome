# load required libraries
library(data.table)

# set working directory
setwd("/xdisk/jneilson/chenyj/base/meta/trait/ags/")

# code is modified from: https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# the length of 35 single-copy genes
cog.length <- read.table(file = "/groups/barberan/metagenome_trait/ags_acn/all_cog_lengths.tsv", header = T,  sep = "\t", row.names = 1)

# read annotation results of 35 single-copy genes
sample.cog.filename <- list.files(pattern = "*single_cogs_count.tsv")

sample.cog <- lapply(sample.cog.filename, function(i){
  i <- paste("", i, sep = "")
  read.table(i,header = T,  sep = "\t", row.names = 1)
  })

# read results of total number of base pairs
sample.totalbp.filename <- list.files(pattern = "*totalbp.txt")

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

# calculate the coverage of the 35 single-copy genes, copy from https://github.com/pereiramemo/AGS-and-ACN-tools/blob/master/code/ags.sh

# in case if coverage of some single-copy genes is zero
for (i in names(sample.cog)) {
  if (sum(!rownames(cog.length) %in% rownames(sample.cog[[i]])) > 0 ) {
    j <- !rownames(cog.length) %in% rownames(sample.cog[[i]])
    missing <- rownames(cog.length)[j]
    cov <- rep(0,length(missing))
    sample.cog[[i]][missing,] <- cov
  }
}

# calculate the coverage of the 35 single-copy genes, that is number of genomes

genome.number <- lapply(sample.cog, function(x) {mean(x[rownames(cog.length), "cov"]/cog.length$value)})

genome.number <- unlist(genome.number)

# calculate ags
ags.result <- sample.totalbp

ags.result$genome_number <- genome.number

ags.result$ags <- (ags.result$totalbp/ags.result$genome_number)/1000000


write.table(ags.result, file = "ags_result.txt", quote=F, sep="\t", col.names = NA)
