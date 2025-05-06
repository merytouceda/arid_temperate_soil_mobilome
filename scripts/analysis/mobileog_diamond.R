# script for preliminary analyses of mobile OG output from BASE

library(tidyverse)
library(vegan)

##############
## BLAST OUTPUT FILTER
##############

mobileog_diamond <- read.table("/Volumes/BunnyBike/mge_urban/base/mobileOG/base_combined_mobileOG__121124.out", header = F)

colnames(mobileog_diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                                "qend", "sstart", "send", "evalue", "bitscore")


# keep only the best hit per gene
mobileog_diamond_filtered <- mobileog_diamond %>%
  dplyr::group_by(qseqid) %>% # group by same read name
  dplyr::arrange(bitscore,.by_group = T) %>% # sort by evalue within group
  slice_head(n = 1) %>% # select the first hit (best one) for each group 
  filter(!pident < 90) %>%
  filter(length > 50)

mean(mobileog_diamond_filtered$length)
summary(mobileog_diamond_filtered$length)
hist(mobileog_diamond_filtered$length)

final_diamond <- mobileog_diamond_filtered  %>%
  separate(sseqid, into = c("mobileOG_id", "gene_name", "other_id", "category", "something", "db", "class"), sep = "\\|") %>%
  mutate(MGE_type = case_when(db %in% c("ISFinder") ~ "Insertion_sequence", 
                              db %in% c("AICE","ICE","CIME","IME","immedb") ~ "Integrative_element",
                              db %in% c("COMPASS","Plasmid") ~ "Plasmid",
                              db %in% c("ACLAME") ~ "Multiple",
                              db %in% c("pVOG","GPD") ~ "Phage")) %>%
  filter(!MGE_type == "Phage") %>%
  filter(!category == "phage") 


plasmid_diamond <- final_diamond %>%
  filter(MGE_type == "Plasmid") 

non_plasmid_diamond <- final_diamond %>%
  filter(!MGE_type == "Plasmid") 

summary(plasmid_diamond$length)
summary(non_plasmid_diamond$length)


# save the contigs
mge_contigs <- non_plasmid_diamond %>%
  select(qseqid)

write.table(mge_contigs, "/Volumes/BunnyBike/mge_urban/base/mobileOG/mges_contigs.txt", quote = F)
