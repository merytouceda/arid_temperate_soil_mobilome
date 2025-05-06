# Genomad output filtering for mge base project

library(tidyverse)
library(vegan)
pal1 = c( "blue", "yellow", "red", "green")


##############
## GENOMAD OUTPUT PARSE & FILTER
##############
setwd("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/all_contigs_plasmid_summary")
files <- list.files(path = "/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/all_contigs_plasmid_summary", pattern="*", full.names=F, recursive=FALSE)

# create empty table to populate
all_samples_genomad <- data.frame(matrix(NA, nrow = 1, ncol = 11))
colnames(all_samples_genomad) <- c("seq_name","length","topology" ,"n_genes","genetic_code" ,"plasmid_score" , "fdr", 
                                   "n_hallmarks","marker_enrichment","conjugation_genes","amr_genes")

for (file in files){
  name <- basename(file)
  # extract sample name from file
  samplename <- str_extract(name, "[^.]*")
  d <- tryCatch(read.table(file, header = T, sep = "\t"), error=function(e) NULL)
  
  # change the name of contig to add the sample name
  d1 <- d %>%
    dplyr::mutate(seq_name = paste(samplename, seq_name, sep = "_"))
  
  # append to form a big dataset
  all_samples_genomad <-rbind(all_samples_genomad, d1)
  
}

# check how many > 10kb
over10kb_summary <- all_samples_genomad%>%
  dplyr::mutate(length_num = as.numeric(length)) %>%
  dplyr::filter(length_num > 10000)


##############
## CIRCULAR BLAST OUTPUT PARSE & FILTER
##############
# open blast files, select those contigs that align to themselves (circular)

# add column names based on blast fmt6
colnames <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
              "qend", "sstart", "send", "evalue", "bitscore")
# do it with all the files: 

# initialize a data frame to add all others and make a huge one
colnames2 <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
               "qend", "sstart", "send", "evalue", "bitscore", "sample_id")
d0 <- data.frame(matrix(NA, nrow = 1, ncol = 13))
colnames(d0) <- colnames2

# list of files: 
setwd("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/blast (1)/")
files <- list.files(path = "/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/blast (1)/", pattern="*_tails.out", full.names=F, recursive=FALSE)

# make a full table of all the samples
# make a full table of all the samples
for (file in files){
  name <- basename(file)
  d <- tryCatch(read.table(file, header = F, sep = "\t"), error=function(e) NULL)
  colnames(d) <- colnames
  # extract sample name from name of file
  samplename <- str_extract(name, "[^.]*")
  
  # create list 
  circular <- d %>%
    filter(qseqid == sseqid) %>%
    mutate(sample_id=samplename) %>%
    mutate(qseqid = paste(samplename, qseqid, sep = "_")) %>%
    mutate(sseqid = paste(samplename, sseqid, sep = "_"))
  
  # append to form a big dataset
  d0 <-rbind(d0, circular)
  
}




##############
## FILTER CIRCULAR OR > 10KB PLASMIDS
##############
# obtain a list of the circular plasmids
circular_plasmids <- d0 %>%
  # filter low quality alignments
  filter(!pident < 95 & !length < 150) %>%
  dplyr::select(qseqid) %>%
  drop_na() %>%
  unique() %>%
  mutate(qseqid = gsub(qseqid, pattern = "_blast_heads_tails_", replacement = "_"))

# filter the output of genomad based on circular plasmids and plasmids over 10kb
all_summary_filtered <- all_samples_genomad %>%   # 61378
  #filter(seq_name %in% circular_plasmids$qseqid) #19
  filter(seq_name %in% circular_plasmids$qseqid | as.numeric(length) > 9999) # 399
  #filter(seq_name %in% circular_plasmids$qseqid & as.numeric(length) > 9999) # 0

#summary(d0$pident)

# write list to filter the contig count table and be able to open it in R
write_csv(all_summary_filtered, "/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/circular_10kb_plasmids.csv")
contig_names <- all_summary_filtered$seq_name
write.table(contig_names, "/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/circular_10kb_plasmids_contignames.csv", sep = "\t", quote = F)

# all sequences genomad says are plasmids
all_contig_names <- all_samples_genomad$seq_name
write.table(all_contig_names, "/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/all_plasmids_contignames.csv", sep = "\t", quote = F)







