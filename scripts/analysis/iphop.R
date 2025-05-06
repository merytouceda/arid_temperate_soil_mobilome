# visualization of iphop results for the base project

library(tidyverse)
library(car)
library(lme4)

pal1 <- c("#D81B60", "#997404", "#004D40")


# open iphop results
iphop <- read_csv("/Volumes/BunnyBike/mge_urban/base/virus/iphop_results/all_host_prediction_to_genus_m90.csv")


# open virus counts and list of filtered virus
virus_counts <- read.table("/Volumes/BunnyBike/mge_urban/base/virus/votus_genomad_count_table.txt", sep = "\t", header = T)
virus_filtered_list <- read.table("/Volumes/BunnyBike/mge_urban/base/virus/filtered_virus_list_phages.txt", sep = " ", header = T)

# make count table of only filtered viruses
virus_count_fixed <- virus_counts %>%
  pivot_longer(-Contig, names_to = "sample", values_to = "count") %>%
  mutate(sample = gsub(sample, pattern = "^X", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  # Change the virus name so the sample is in front and then the contig (like in virus_filtered_list)
  mutate(virus_sample = str_extract(Contig, pattern = "_[0-9]*$")) %>%
  mutate(virus_sample = gsub(virus_sample, pattern = "_", replacement = "")) %>%
  mutate(virus_contig = str_extract(Contig, pattern = "k[0-9]*_[0-9]*_")) %>%
  mutate(virus_contig = gsub(virus_contig, pattern = "_$", replacement = "")) %>%
  mutate(virus_name = paste(virus_sample, virus_contig , sep = "_")) %>%
  dplyr::select(-c("Contig", "virus_sample", "virus_contig")) %>%
  # filter out those viruses that are not in the filtered list
  filter(virus_name %in% virus_filtered_list$virus_name) %>%
  pivot_wider(names_from = "sample", values_from = "count")

# load metadata and match with counts
md <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_votus.csv")


virus_counts_long <- virus_count_fixed  %>%
  pivot_longer(-virus_name, names_to = "Sample", values_to = "Counts") 

md_filtered <- md %>%
  #mutate(Sample = gsub(sample_id, pattern = "_[0-9]", replacement = "")) %>%
  filter(Sample %in% virus_counts_long$Sample) %>%
  filter(!general_env_feature == "Australia") %>%
  filter(!is.na(vegetation_type)) %>%
  filter(!vegetation_type == "Forest")

virus_counts_filtered <- virus_counts_long %>%
  filter(Sample %in% md_filtered$Sample) %>%
  pivot_wider(names_from = virus_name, values_from = Counts ) %>%
  arrange(desc(as.numeric(Sample))) %>%
  column_to_rownames(var = "Sample") 


md_final<- md_filtered %>%
  arrange(desc(Sample)) %>%   
  column_to_rownames(var = "Sample") 

# make sure the rows are in the same order
rownames(virus_counts_filtered) == rownames(md_final)


# --------------------------------------------------------------------------------------------------------------------
##############
# RPKM
##############
# open length file
virus_filtered_len <- read.table("/Volumes/BunnyBike/mge_urban/base/virus/filtered_virus_length_phages.txt", header = T)

# select the total reads per sample from the metadata
total_reads <- md_final %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(Sample, tot_reads)

# 2. pivot longer the count table to make every gene in every sample a row
virus_long <- virus_counts_filtered %>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = -c("Sample"), values_to = "Count", names_to= "virus_name") 

# 3. join the total reads to each row based on sample
virus_long_plusreads <- left_join(virus_long, total_reads, by = "Sample")
# 4. join the gene length to each row based on gene
virus_long_plusreads_plustlen <- left_join(virus_long_plusreads, virus_filtered_len, by = "virus_name")

# 5. compute the RPKM calculation and create count table from it 
rpkm_virus_final <- virus_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Count*10e6)/(tot_reads*as.numeric(length))) %>%
  dplyr::select(-c("Count", "tot_reads", "length")) %>%
  #arrange(Sample) %>% # run with this option only for Mantel tests! 
  pivot_wider(names_from = Sample, values_from = RPKM)%>%
  column_to_rownames(var = "virus_name") %>%
  unique() %>%
  drop_na()

# ------------------------------------------------------------------------------------------------
# prepare iphop output

iphop_clean <- iphop %>%
  # filter problematic rows and header rows from cat
  filter(!str_detect(Virus, "provirus")) %>%
  filter(!Virus == "Virus") %>%
  # first we fix the virus name so it matches the count table format
  separate(Virus, into = c("contig1", "contig2", "sample"), sep = "_") %>%
  mutate(virus_name = paste(sample, contig1, contig2, sep = "_")) %>%
  dplyr::select(-c("contig1", "contig2", "sample")) %>%
  # now we deal with the taxonomy
  separate(`Host genus`, into = c("domain", "phylum", "class", "order", "family", "genus"), sep = ";")


iphop_join <- iphop_clean %>%
  select(c("virus_name", "phylum")) %>%
  distinct(virus_name, .keep_all = T)


# ------------------------------------------------------------------------------------------------
## Join with iphop and count the hosts

md_final_join <- md_final %>%
  rownames_to_column(var = "Sample")

counts_pab_iphop <- rpkm_virus_final %>%
  ungroup() %>%
  rownames_to_column(var = "virus_name") %>%
  pivot_longer(!virus_name, names_to = "Sample", values_to = "pab") %>%
  left_join(iphop_join, by = "virus_name") %>%
  drop_na() %>%
  left_join(., md_final_join, by = "Sample") %>%
  group_by(phylum, general_env_feature, vegetation_type) %>%
  summarize(counts = sum(pab)) %>%
  mutate(pab = case_when(counts == 0 ~ 0, 
                         TRUE ~ 1)) %>%
  mutate(combined_var = paste(general_env_feature, vegetation_type)) %>%
  arrange(desc(counts)) %>%
  # calculate percentage to plot
  group_by(combined_var) %>%
  mutate(percentage = (counts / sum(counts)) * 100) %>%
  ungroup()%>%
  arrange(combined_var, desc(percentage))
  

  
pal3 <- c("orange" ,"#FEE08B", "#D53E4F" ,"#3288BD")
  
ggplot(counts_pab_iphop, aes(x = vegetation_type, y = percentage, fill = phylum, group = counts))+
  geom_col()+
  scale_fill_manual(values = pal3)+
  xlab("")+
  ylab("Proportion of counts (%)")+
  facet_wrap(~general_env_feature)+
  theme_bw()+
  theme(legend.position = "none")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/virus_hosts_phyla_nolegend.pdf", device = "pdf", width = 6, height = 5 , units = "in")
  
# have the taxonomy of 45 out of 337 virus (13%)




