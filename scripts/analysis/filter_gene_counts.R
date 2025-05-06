# script to filter out the gene count table for base

library(tidyverse)


# load data
## Input files and directories
args<-commandArgs(TRUE)

gene_count_table=args[1]

gene_count_table <- read.table("Downloads/head_all_genes_count_table.txt", sep = "\t", header = T)

# Change low coverage read counts for 0
# Count = 0 if coverage < 100X
gene_length 

#counts*150/gene_length = coverage
#min coverage = 100X



#Now filter by prevalence
gene_count_table_prvalence <- gene_count_table %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "count") %>%
  mutate(pab = case_when(count < 0 ~ 0,  # if count < 50 then coverage is too low
                         TRUE ~ 1)) %>%
  group_by(Contig) %>%
  summarize(pab_total = sum(pab)) %>%
  filter(pab_total > 19). # change this for argument n*0.1

#filter the table
filter_gene_count_table <- gene_count_table %>%
  filter(Contig %in% gene_count_table_prvalence$Contig)


# Abundance filter (optional)
# calculate summary of 
filter_gene_count_table_long <- filter_gene_count_table%>%
  pivot_longer(-Contig, names_to = "sample", values_to = "counts") %>%
  group_by(Contig) %>%
  summarize(suma = sum(counts)) %>%
  filter(suma > 100)

filter_gene_count_table_final <- filter_gene_count_table %>%
  filter(Contig %in% filter_gene_count_table_long$Contig)



write_table(filter_gene_count_table, "filtered_all_gene_counts.txt", quote = F)
