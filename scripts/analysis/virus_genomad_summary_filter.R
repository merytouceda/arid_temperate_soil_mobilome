# virus genomad summary visualizations for selecting viruses

# According to Graham et.al. 2024: 

# length > 1k
# (1) contigs of at least 1 kb with high similarity to genomes in the CheckV database 
#(that is, that had high- or medium-quality completeness estimates) or that contained direct terminal repeats 
#were automatically selected; 

# 10k > length > 5k
# (2) contigs shorter than 10 kb and longer than 5 kb were required to have a geNomad virus score higher than 0.9, 
#to encode at least one virus hallmark AND to have a virus marker enrichment higher than 2.0. 

# length > 10k
# (3) contigs longer than 10 kb were required to have a geNomad virus score higher than 0.8 
#and to EITHER encode one virus hallmark (for example, terminase, capsid proteins, portal protein and so on), 
#as determined by geNomad,OR to have a geNomad virus marker of at least 5.0; 

# 337 vOTUs of bacteriophages

library(tidyverse)
library(ggplot2)
library(see)
library(vegan)
library(lme4)
library(car)


pal1 <- c("#D81B60", "#FFC107", "#004D40")
pal2 <- c("#D81B60", "#997404", "#004D40")


colnames(virus_summary)

'''
############################
### VIRUS SUMMARY 
############################

# Open the virus summary files one by one and make a table
path <- "/Volumes/BunnyBike/mge_urban/base/virus/all_contigs_virus_summary"
setwd("/Volumes/BunnyBike/mge_urban/base/virus/all_contigs_virus_summary")

# initialize empty d0
columns <- c("seq_name","length", "topology" ,"coordinates","n_genes","genetic_code","virus_score",
             "fdr","n_hallmarks", "marker_enrichment","taxonomy", "sample_id", "virus_name")

d0 <- data.frame(matrix(NA, nrow = 1, ncol = 13))
colnames(d0) <- columns

# list of files: 
files <- list.files(path = path, pattern="*.tsv", full.names=F, recursive=FALSE)

ok # make a full table of all the samples
for (file in files){
  
  # extract the sample number to add to the virus name later on
  name <- basename(file)
  samplename <- str_extract(name, "[^\\.]*")
  
  # open the file
  d <- tryCatch(read.table(file, header = T, sep = "\t"), error=function(e) NULL)



  d <- d %>%
    mutate(sample_id=samplename) %>%
    mutate(virus_name = paste(sample_id, seq_name, sep = "_"))
  
  # append to form a big dataset
  d0 <-rbind(d0, d)
  
}


# make virus_summary table from the d0 table
virus_summary <- d0 %>%
  filter(!is.na(seq_name)) %>%
  dplyr::select(-c("seq_name", "sample_id")) %>%
  filter(!str_detect(virus_name, "all_virus")) %>%
  # keep only the ones we know are bacteriophages
  filter(str_detect(taxonomy, "Uroviricota") | str_detect(taxonomy, "Phixviricota") | str_detect(taxonomy, "Hofneiviricota") )


###############################
### GENES
###############################

path <- "/Volumes/BunnyBike/mge_urban/base/virus/all_contigs_virus_genes_summary"
setwd("/Volumes/BunnyBike/mge_urban/base/virus/all_contigs_virus_genes_summary")

# initialize empty d0
columns <- c("gene", "start", "end","length","strand","gc_content","genetic_code", "rbs_motif","marker" ,"evalue","bitscore",
             "uscg" , "plasmid_hallmark", "virus_hallmark" ,"taxid", "taxname" , "annotation_conjscan" , "annotation_amr",
             "annotation_accessions","annotation_description", "sample_id", "gene_name")


d0 <- data.frame(matrix(NA, nrow = 1, ncol = 22))
colnames(d0) <- columns

# list of files: 
files <- list.files(path = path, pattern="*.tsv", full.names=F, recursive=FALSE)

# make a full table of all the samples
for (file in files){
  
  # extract the sample number to add to the virus name later on
  name <- basename(file)
  samplename <- str_extract(name, "[^\\.]*")
  
  # open the file
  d <- tryCatch(read.table(file, header = T, sep = "\t"), error=function(e) NULL)
  
  
  d <- d %>%
    mutate(sample_id=samplename) %>%
    mutate(gene_name = paste(sample_id, gene, sep = "_"))
  
  # append to form a big dataset
  d0 <-rbind(d0, d)
  
}



genes_summary <- d0 %>%
  filter(!is.na(gene_name)) %>%
  dplyr::mutate(virus_name = gsub(gene_name, pattern = "_[0-9]*$", replacement = "")) %>%
  dplyr::select(-c("sample_id"))
  

gene_summary_tojoin <- genes_summary %>%
  dplyr::select(c("virus_name", "virus_hallmark")) %>% 
  dplyr::group_by(virus_name) %>%
  dplyr::summarize(hallmark = sum(virus_hallmark))


virus_summary_plusgene <- virus_summary %>%
  left_join(., gene_summary_tojoin, by = "virus_name")

# the hallmark genes (n_hallmarks) in the summary virus file is already the result of this, so it is not necessary to do it


###############################
### CHECKV
###############################

path <- "/Volumes/BunnyBike/mge_urban/base/virus/all_checkv_summary"
setwd("/Volumes/BunnyBike/mge_urban/base/virus/all_checkv_summary")

# initialize empty d0
columns <- c("contig_id","contig_length","provirus" ,"proviral_length" , "gene_count" ,"viral_genes",
             "host_genes", "checkv_quality","miuvig_quality","completeness" ,"completeness_method" ,"contamination" ,
             "kmer_freq","warnings" ,"sample_id", "virus_name")


d0 <- data.frame(matrix(NA, nrow = 1, ncol = 16))
colnames(d0) <- columns

# list of files: 
files <- list.files(path = path, pattern="*.tsv", full.names=F, recursive=FALSE)

# make a full table of all the samples
for (file in files){
  
  # extract the sample number to add to the virus name later on
  name <- basename(file)
  samplename <- str_extract(name, "[^_]*")
  
  # open the file
  d <- tryCatch(read.table(file, header = T, sep = "\t"), error=function(e) NULL)
  
  
  d <- d %>%
    mutate(sample_id=samplename) %>%
    mutate(virus_name = paste(sample_id, contig_id, sep = "_"))
  
  # append to form a big dataset
  d0 <-rbind(d0, d)
  
}


checkv_summary <- d0 %>%
  filter(!is.na(virus_name)) %>%
  dplyr::select(-c("sample_id"))


checkv_summary_tojoin <- checkv_summary %>%
  dplyr::select(c("virus_name", "checkv_quality"))


virus_summary_checkv<- virus_summary %>%
  left_join(., checkv_summary_tojoin, by = "virus_name")

# the hallmark genes (n_hallmarks) in the summary virus file is already the result of this, so it is not necessary to do it




###########
### FILTERING
###########
# see length distributions
summary(as.numeric(virus_summary$length))
ggplot(virus_summary_1kb, aes(x = as.numeric(length)))+
  geom_density()



  
# (1) 
# length > 1k
# (1) contigs of at least 1 kb with high similarity to genomes in the CheckV database 
#(that is, that had high- or medium-quality completeness estimates) or that contained direct terminal repeats 
virus_summary_1kb_graham <- virus_summary_checkv %>%
  filter(as.numeric(length) > 1000)  %>%
  filter(as.numeric(length) < 5000) %>%
  filter(topology %in% c("DTR" ,"ITR", "Provirus") |
           checkv_quality %in% c("Complete","High-quality", "Medium-quality")) # high similarity to genomes in the CheckV database or that contained direct terminal repeats 

virus_summary_1kb_graham_names <- virus_summary_1kb_graham %>%
  dplyr::select(virus_name)

# (2) 
# 10k > length > 5k
# (2) contigs shorter than 10 kb and longer than 5 kb were required to have a geNomad virus score higher than 0.9, 
#to encode at least one virus hallmark AND to have a virus marker enrichment higher than 2.0. 

virus_summary_5kb <- virus_summary %>%
  filter(as.numeric(length) > 5000)

virus_summary_5kb_graham <- virus_summary %>%
  filter(as.numeric(length) > 5000)%>%
  filter(as.numeric(length) < 10000)%>%
  filter(as.numeric(virus_score) > 0.9) %>%  # a geNomad virus score higher than 0.9
  filter(as.numeric(marker_enrichment) > 2.0 &&  n_hallmarks > 0)  # encode at least one virus hallmark AND to have a virus marker enrichment higher than 2.0

virus_summary_5kb_graham_names <- virus_summary_5kb_graham %>%
  dplyr:: select(virus_name)


# (3) 
# length > 10k
# (3) contigs longer than 10 kb were required to have a geNomad virus score higher than 0.8 
#and to EITHER encode one virus hallmark (for example, terminase, capsid proteins, portal protein and so on), 
#as determined by geNomad,OR to have a geNomad virus marker of at least 5.0; 
virus_summary_10kb <- virus_summary %>%
  filter(as.numeric(length) > 10000)

virus_summary_10kb_graham <- virus_summary %>%
  filter(as.numeric(length) > 10000) %>%
  filter(as.numeric(virus_score) > 0.8) %>% # geNomad virus score higher than 0.8 
  filter(as.numeric(marker_enrichment) > 5.0 | n_hallmarks > 0 ) # to EITHER encode one virus hallmark OR to have a geNomad virus marker of at least 5.0

virus_summary_10kb_graham_names <- virus_summary_10kb_graham %>%
  dplyr::select(virus_name)


virus_filtered_list <- rbind(virus_summary_10kb_graham_names , virus_summary_5kb_graham_names, virus_summary_1kb_graham_names)

virus_filtered_list <-virus_filtered_list  %>%
  dplyr::mutate(virus_name = gsub(virus_name, pattern = "|provirus_[0-9]*_[0-9]*", replacement = "")) %>%
  dplyr::mutate(virus_name = gsub(virus_name, pattern = "\\|$", replacement = ""))

# Virus list total = 418 
write.table(virus_filtered_list, "/Volumes/BunnyBike/mge_urban/base/virus/filtered_virus_list_phages.txt", quote = F)


# get viral length
virus_summary_10kb_len<- virus_summary_10kb_graham %>%
  dplyr::select(c("virus_name", "length"))

virus_summary_1kb_len<- virus_summary_1kb_graham %>%
  dplyr::select(c("virus_name", "length"))

virus_summary_5kb_len<- virus_summary_5kb_graham %>%
  dplyr::select(c("virus_name", "length"))

virus_filtered_len <- rbind(virus_summary_10kb_len , virus_summary_5kb_len, virus_summary_1kb_len)

virus_filtered_len <- virus_filtered_len %>%
  mutate(virus_name = gsub(virus_name, pattern = "|provirus_[0-9]*_[0-9]*", replacement = "")) %>%
  mutate(virus_name = gsub(virus_name, pattern = "\\|$", replacement = ""))

write.table(virus_filtered_len, "/Volumes/BunnyBike/mge_urban/base/virus/filtered_virus_length_phages.txt", quote = F)

'''

#### #### #### #### #### #### #### #### #### #### #### 
#### COUNTS
#### #### #### #### #### #### #### #### #### #### #### 

#virus_counts <- read.table("/Volumes/BunnyBike/mge_urban/base/virus/virus_genomad_count_table.txt", sep = "\t", header = T)
virus_counts <- read.table("/Volumes/BunnyBike/mge_urban/base/virus/votus_genomad_count_table.txt", sep = "\t", header = T)
#virus_filtered_list <- read.table("/Volumes/BunnyBike/mge_urban/base/virus/filtered_virus_list.txt", sep = " ", header = T)
virus_filtered_list <- read.table("/Volumes/BunnyBike/mge_urban/base/virus/filtered_virus_list_phages.txt", sep = " ", header = T)


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

# load metadata
#md <- read.table("/Volumes/BunnyBike/mge_urban/base/base_metadata.txt", sep = "\t", header = T)
md <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_ptus.csv")

virus_counts_long <- virus_count_fixed  %>%
  pivot_longer(-virus_name, names_to = "Sample", values_to = "Counts") 

md_filtered <- md %>%
  #mutate(Sample = gsub(sample_id, pattern = "_[0-9]", replacement = "")) %>%
  filter(Sample %in% virus_counts_long$Sample) %>%
  filter(!general_env_feature == "Australia") %>%
  filter(!is.na(vegetation_type)) %>%
  filter(!vegetation_type == "Forest") %>%
  filter(!Sample == "8459")

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
### RICHNESS
summary(rowSums(virus_counts_filtered))
hist(rowSums(virus_counts_filtered))

# RICHNESS
md_final$virus_rich <- specnumber(rrarefy(virus_counts_filtered, sample =  700))

# SHANNON with rarefaction
md_final$virus_shannon <- diversity(rrarefy(virus_counts_filtered, sample =  700), index = "shannon")

'''
# SHANNON with normalization
shannon_diversity <- function(x) {
x <- x[x > 0]  # Remove zeros
p <- x/sum(x)
-sum(p * log(p))
}

shannon_values <- as.data.frame(apply(virus_counts_filtered, 1, shannon_diversity)) %>%
  rownames_to_column(var = "Sample")
colnames(shannon_values) <- c("Sample", "shannon_norm")

md_final <- md_final %>%
  left_join(., shannon_values, by = "Sample")

# Convert to effective number of species (Hill number, q=1)
exp_shannon <- exp(shannon_values)


# SHANNON with small sample bias correction
shannon_corrected <- function(x) {
  S <- sum(x > 0)  # Number of species observed
  n <- sum(x)      # Total count
  H <- shannon_diversity(x)
  H + (S-1)/(2*n)  # Bias correction
}

shannon_corr_values <- as.data.frame(apply(virus_counts_filtered, 1, shannon_corrected)) %>%
  rownames_to_column(var = "Sample")
colnames(shannon_corr_values) <- c("Sample", "shannon_samplebias")

md_final <- md_final %>%
  left_join(., shannon_corr_values, by = "Sample")

# simpson
md_final$virus_simpson <- diversity(rrarefy(virus_counts_filtered, sample =  700), index = "simpson")

'''

ggplot(md_final, aes(x = vegetation_type, y = virus_rich))+
  #geom_violin(aes(fill =urban.natural, color = urban.natural), binwidth = 1, alpha=0.4)+
  geom_boxplot(alpha=0.6, outlier.shape = NA )+
  geom_jitter(aes(color= vegetation_type), size = 2.5)+
  xlab(NULL) + 
  ylab("Virus richness")+
  scale_fill_see()+
  scale_color_see()+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~general_env_feature)
ggsave("/Volumes/BunnyBike/mge/base/figures/virus_rich.pdf", device = "pdf", width = 5, height = 5 , units = "in")

# nice plot
ggplot(md_final, aes(x = vegetation_type, y = virus_rich, color = vegetation_type, fill = vegetation_type))+
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA, linewidth = 0.5) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 2)+
  scale_fill_manual(values = pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Virus richness")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/votus_rich_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")

Anova(lmer(virus_rich ~  vegetation_type * general_env_feature + (1|season), data = md_final), type = 3)
#vegetation_type                      0.9264  2    0.62928    
#general_env_feature                  3.2956  1    0.06946 .  
#vegetation_type:general_env_feature  6.7430  2    0.03434 * 
Anova(lmer(virus_shannon ~  vegetation_type * general_env_feature + (1|season), data = md_final), type = 3)
#vegetation_type                       0.9027  2     0.6368    
#general_env_feature                   0.0739  1     0.7858    
#vegetation_type:general_env_feature   4.5986  2     0.1003
r2(lmer(virus_rich~  vegetation_type * general_env_feature + (1|season), data = md_final))
r2(lmer(virus_shannon~  vegetation_type * general_env_feature + (1|season), data = md_final))

# alternative aesthetic

md_final$combined_var <- paste(md_final$general_env_feature, md_final$vegetation_type)
md_final$combined_var <- factor(md_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                              "Arid Shrubland", "Temperate Shrubland", 
                                                                              "Arid Woodland", "Temperate Woodland"))

ggplot(md_final, aes(x = vegetation_type, y = virus_rich, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  ylab("Virus richness")+ 
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/votus_rich_paralell.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# Plot to check if diversity correlates with sample size
# Create a data frame with your diversity metrics and sample sizes
plot_data <- data.frame(
  SampleSize = rowSums(virus_counts_filtered),
  SimpsonDiv = md_final$virus_rich
)

# Create the plot
ggplot(plot_data, aes(x = SampleSize, y = SimpsonDiv)) +
  geom_point() +                                           # Add points
  geom_smooth(method = "lm", se = TRUE, color = "red") +  # Add regression line with confidence interval
  labs(
    x = "Sample size (number of counts)",
    y = "Richness",
  ) +
  theme_bw()
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/rich_samplesize.pdf", device = "pdf", width = 5, height = 5 , units = "in")







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

# --------------------------------------------------------------------------------------------------------------------
############################
# beta-diversity
############################
virus.bray <- vegdist(t(rpkm_virus_final), method="bray")

#ordination (non-multidimensional scaling)
virus.nmds <- metaMDS(virus.bray, k=2, try = 100)
md_final$Axis01 = virus.nmds$points[,1]
md_final$Axis02 = virus.nmds$points[,2]
virus.nmds$stress 

ggplot(md_final, aes(Axis01, Axis02))+
  geom_point(aes(color=vegetation_type, shape = general_env_feature), size=4)+
  #stat_ellipse(aes(color = vegetation_type))+
  stat_ellipse(aes(group = general_env_feature, linetype = general_env_feature))+
  scale_shape_manual(values=c(16, 13))+
  scale_color_manual(values = pal1)+
  theme_bw()+
  theme(legend.position="none", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/votus_beta_rpkm.pdf", device = "pdf", width = 7, height = 5 , units = "in")

adonis2(virus.bray ~ season + vegetation_type * general_env_feature, data = md_final,  permutations = 999, method = "bray")
#vegetation_type                       2    1.557 0.03314 2.7026  0.001 ***
#  general_env_feature                   1    2.184 0.04648 7.5807  0.001 ***
#  vegetation_type:general_env_feature   2    1.186 0.02525 2.0590  0.003 *


# soil properties
adonis2(virus.bray ~ vegetation_type + general_env_feature + ph, data = md_final,  permutations = 999, method = "bray")
#ph                    1    0.729 0.01550 2.5187  0.002 ** 
adonis2(virus.bray ~ vegetation_type + general_env_feature + organic_carbon, data = md_final,  permutations = 999, method = "bray")
#organic_carbon        1    1.220 0.02596 4.2663  0.001 **


md_final$combined_var <- paste(md_final$general_env_feature, md_final$vegetation_type)
md_final$combined_var <- factor(md_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                  "Arid Shrubland", "Temperate Shrubland", 
                                                                  "Arid Woodland", "Temperate Woodland"))
# new aesthetics: 
ggplot(md_final, aes(Axis01, Axis02))+
  #geom_point(aes(color=vegetation_type, shape = general_env_feature), size=4)+
  geom_point(aes(color=combined_var, shape = general_env_feature), size=4)+
  #stat_ellipse(aes(color = vegetation_type))+
  stat_ellipse(aes(group = general_env_feature, linetype = general_env_feature))+
  #scale_shape_manual(values=c(16, 13))+
  scale_shape_manual(values=c(18, 16))+
  #scale_color_manual(values = pal2)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  theme_bw()+
  theme(legend.position="none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/votus_beta_rpkm_combinedvar.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# bray curtis dissimilarity separate by arid and temperate
md_arid <- md_final %>%
  filter(general_env_feature == "Arid") %>%
  mutate(sample = gsub(sample_id, pattern = "_[0-9]$", replacement = ""))

rpkm_virus_final_arid <- rpkm_virus_final %>% 
  select(md_arid$sample)

rpkm_virus_final_temperate <- rpkm_virus_final %>% 
  select(!md_arid$sample)

virus.bray.arid <- vegdist(t(rpkm_virus_final_arid), method="bray")
virus.bray.temperate <- vegdist(t(rpkm_virus_final_temperate), method="bray")


# --------------------------------------------------------------------------------------------------------------------
############################
# abundance
############################
sum_rel_ab_virus <- as.data.frame(t(rpkm_virus_final)) 

sum_rel_ab_virus_final <- sum_rel_ab_virus %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(virus_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(virus_relab) %>%
  rownames_to_column(var = "Sample")

md_final<- md_final %>%
  #rownames_to_column(var = "Sample") %>%
  left_join(., sum_rel_ab_virus_final, by = "Sample")

ggplot(md_final, aes(x = vegetation_type, y = virus_relab))+
  #geom_violin(aes(fill =urban.natural, color = urban.natural), binwidth = 1, alpha=0.4)+
  geom_boxplot(alpha=0.6, outlier.shape = NA )+
  geom_jitter(aes(color= vegetation_type), size = 2.5)+
  xlab(NULL) + 
  ylab("Plasmid relative abundance")+
  scale_fill_manual(values = pal1)+
  scale_color_manual(values = pal1)+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~general_env_feature)
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/virus_relab_veg.pdf", device = "pdf", width = 5, height = 5 , units = "in")


ggplot(md_final, aes(x = vegetation_type, y = virus_relab, color = vegetation_type, fill = vegetation_type))+
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA, linewidth = 0.5) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 2)+
  scale_fill_manual(values = pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Virus relative abundance (RPKM)")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/votus_relab_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")


anova(lm(virus_relab ~  vegetation_type * general_env_feature, data = md_final))
Anova(lmer(virus_relab ~  vegetation_type * general_env_feature + (1|season), data = md_final), type = 3)
#vegetation_type                     0.1774  2   0.915131   
#general_env_feature                 1.9610  1   0.161407   
#vegetation_type:general_env_feature 2.8464  2   0.240939 



# --------------------------------------------------------------------------------------------------------------------
## Length
counts_pab <- rpkm_virus_final
counts_pab[counts_pab > 0] <- 1

md_final_join <- md_final %>%
  rownames_to_column(var = "Sample")

counts_pab_length <- counts_pab %>%
  ungroup() %>%
  rownames_to_column(var = "virus_name") %>%
  pivot_longer(!virus_name, names_to = "Sample", values_to = "pab") %>%
  left_join(virus_filtered_len, by = "virus_name") %>%
  filter(pab == 1) %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarize(virus_av_len = mean(length, na.rm = TRUE), 
                   virus_min_len = min(length, na.rm = TRUE), 
                   virus_max_len = max(length, na.rm = TRUE)) %>%
  left_join(md_final_join, by = "Sample")


ggplot(counts_pab_length, aes(x = vegetation_type, y = virus_av_len, color = vegetation_type, fill = vegetation_type))+
  geom_boxplot(
    width = .2, fill = "white",
    size = 1.5, outlier.shape = NA, linewidth = 0.5) +
  ggdist::stat_halfeye(
    adjust = .33, ## bandwidth
    width = .67, 
    color = NA, ## remove slab interval
    position = position_nudge(x = .15)) +
  gghalves::geom_half_point(
    side = "l", 
    range_scale = .3, 
    alpha = .5, size = 2)+
  scale_fill_manual(values =pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Virus length")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/virus_length_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")

Anova(lmer(virus_av_len ~  vegetation_type * general_env_feature  + (1|season), data = counts_pab_length), type = 3)
#vegetation_type                       3.4390  2     0.1792    
#general_env_feature                   2.2592  1     0.1328    
#vegetation_type:general_env_feature   0.7000  2     0.7047
r2(lmer(virus_av_len~  vegetation_type * general_env_feature + (1|season), data = counts_pab_length))

# nice plot
counts_pab_length$combined_var <- paste(counts_pab_length$general_env_feature, counts_pab_length$vegetation_type)
counts_pab_length$combined_var <- factor(counts_pab_length$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                              "Arid Shrubland", "Temperate Shrubland", 
                                                              "Arid Woodland", "Temperate Woodland"))
ggplot(counts_pab_length, aes(x = vegetation_type, y = virus_av_len, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  xlab(NULL) + 
  ylab("Average plasmid size")+
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/votus_length_nice_parallel.pdf", device = "pdf", width = 5, height = 5 , units = "in")




# ------------------------------------------------------------------------------------------------------------------------------
# length frequency distributions

counts_pab <- rpkm_virus_final
counts_pab[counts_pab > 0] <- 1

md_final_final_join <- md_final %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(c("Sample", "general_env_feature", "vegetation_type"))

counts_pab_length <- counts_pab %>%
  ungroup() %>%
  rownames_to_column(var = "virus_name") %>%
  pivot_longer(!virus_name, names_to = "Sample", values_to = "pab") %>%
  left_join(virus_filtered_len, by = "virus_name") %>%
  filter(pab == 1) %>%
  left_join(md_final_final_join, by = "Sample")


ggplot(counts_pab_length, aes(x = length, fill = vegetation_type))+
  geom_histogram(position = "identity", alpha = 0.6, bins = 50) +
  facet_wrap(~general_env_feature)+
  scale_fill_manual(values = pal2)+
  ylab("Frequency") + 
  xlab("Virus length")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/votus_length_distribution.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# --------------------------------------------------------------------------------------------------------------------
# Viral lifestyle
phatyp <- read_tsv("/Volumes/BunnyBike/mge_urban/base/virus/phatyp_prediction.tsv")

phatyp_clean <- phatyp %>%
  distinct(Accession, .keep_all = T) %>%
  filter(!str_detect(Accession, "provirus")) %>% # there are only 6 provirus, I rather take them out than deal with them
  separate(Accession, into = c("Contig1", "Contig2", "Sample"), sep = "_") %>%
  unite("virus_name", c("Sample", "Contig1", "Contig2"), sep = "_")

### NOTE: fix the virus_names so they match the counts virus_names

counts_pab <- rpkm_virus_final
counts_pab[counts_pab > 0] <- 1

md_final_final_join <- md_final %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(c("Sample", "general_env_feature", "vegetation_type", "season"))


counts_pab_phatyp <- counts_pab %>%
  ungroup() %>%
  rownames_to_column(var = "virus_name") %>%
  pivot_longer(!virus_name, names_to = "Sample", values_to = "pab") %>%
  left_join(phatyp_clean, by = "virus_name") %>%
  #filter(pab == 1) %>%
  left_join(md_final_final_join, by = "Sample") %>%
  filter(TYPE %in% c("virulent", "temperate")) 

# how many viruses were annotated
counts_pab_phatyp_howmany <- counts_pab_phatyp %>%
  count(virus_name)
# 265/337 viruses annotated (0.79)



# calculate ratio
counts_pab_phatyp_ratio <- counts_pab_phatyp %>%
  # first count how many viruses of each type in each sample
  group_by(Sample, TYPE) %>%
  summarize(counts = sum(pab)) %>%
  # make types two columns and replace NAs for 0s
  pivot_wider(names_from = "TYPE", values_from = "counts") %>%
  mutate(across(everything(), ~replace_na(., 0))) %>%
  mutate(ratio_vt = virulent/(temperate+virulent)) %>%
  left_join(md_final_final_join, by = "Sample")

# percentage of virulent
mean(counts_pab_phatyp_ratio$ratio_vt)
sd(counts_pab_phatyp_ratio$ratio_vt)

counts_pab_phatyp_ratio_arid <- counts_pab_phatyp_ratio %>%
  filter(general_env_feature == "Arid")
mean(counts_pab_phatyp_ratio_arid$ratio_vt)

counts_pab_phatyp_ratio_temperate <- counts_pab_phatyp_ratio %>%
  filter(!general_env_feature == "Arid")
mean(counts_pab_phatyp_ratio_temperate$ratio_vt)



# Plot
counts_pab_phatyp_ratio$combined_var <- paste(counts_pab_phatyp_ratio$general_env_feature, counts_pab_phatyp_ratio$vegetation_type)
counts_pab_phatyp_ratio$combined_var <- factor(counts_pab_phatyp_ratio$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                                    "Arid Shrubland", "Temperate Shrubland", 
                                                                                    "Arid Woodland", "Temperate Woodland"))
ggplot(counts_pab_phatyp_ratio, aes(x = vegetation_type, y = ratio_vt, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  xlab(NULL) + 
  ylab("Proportion of virulent viruses")+
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/virulent_ratio_number.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(ratio_vt ~  vegetation_type * general_env_feature  + (1|season), data = counts_pab_phatyp_ratio), type = 3)
#vegetation_type                       1.7943  2     0.4077    
#general_env_feature                   2.3383  1     0.1262    
#vegetation_type:general_env_feature   0.3727  2     0.8300
r2(lmer(ratio_vt~  vegetation_type * general_env_feature + (1|season), data = counts_pab_phatyp_ratio))




# --------------------------------------------------------------------------------------------------------------------
## Tukey test
library(emmeans)
library(multcomp)

m1 <- lmer(virus_shannon ~  vegetation_type * general_env_feature  + (1|season), data = md_final)
m3 <- lmer(virus_av_len ~  vegetation_type * general_env_feature  + (1|season), data = counts_pab_length)

# Step 1: Get estimated marginal means for the interaction
emm_interaction <- emmeans(m1, ~ vegetation_type * general_env_feature)

# Step 2: Perform pairwise comparisons with Tukey adjustment
pairs_result <- pairs(emm_interaction, adjust = "tukey")

# Step 3: View the results
summary(pairs_result)


# Optional: Create a compact letter display
cld_result <- cld(emm_interaction, Letters = letters, adjust = "tukey")
cld_result


emm_df <- as.data.frame(emm_interaction)

# Create the interaction plot
ggplot(emm_df, aes(x = general_env_feature, y = emmean, color = vegetation_type, group = vegetation_type)) +
  geom_point(size = 3) +
  geom_line(linewidth = 1) +
  geom_errorbar(aes(ymin = emmean - SE, ymax = emmean + SE), width = 0.1) +
  labs(x = "Climatic Classification", 
       y = "Virus richness") +
  theme_bw() +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values = pal1)
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/virus_rich_interaction.pdf", device = "pdf", width = 3, height = 5 , units = "in")


# --------------------------------------------------------------------------------------------------------------------
# save md
av_len <- counts_pab_length %>%
  select(c("Sample", "virus_av_len"))

counts_pab_phatyp_join <- counts_pab_phatyp_ratio %>%
  select(c("Sample", "ratio_vt"))

md_final_final <- md_final %>%
  rownames_to_column(var = "Sample") %>%
  left_join(.,av_len, by = "Sample") %>%
  left_join(., counts_pab_phatyp_join, by = "Sample")

write_csv(md_final_final, "/Volumes/BunnyBike/mge_urban/base/md_final_final_votus.csv")



############################
# beta-diversity separated arid and temperate
############################
md_arid <- md_final_final %>%
  rownames_to_column(var = "Sample") %>%
  filter(general_env_feature == "Arid")

rpkm_virus_arid <- rpkm_virus_final %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "Count") %>%
  filter(Sample %in% md_arid$Sample) %>%
  pivot_wider(names_from = "Sample", values_from = "Count") %>%
  column_to_rownames(var = "Contig")

rpkm_virus_temperate <- rpkm_virus_final %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "Count") %>%
  filter(!Sample %in% md_arid$Sample) %>%
  pivot_wider(names_from = "Sample", values_from = "Count") %>%
  column_to_rownames(var = "Contig")

virus.bray_arid <- vegdist(t(rpkm_virus_arid), method="bray")
virus.bray_temperate <- vegdist(t(rpkm_virus_temperate), method="bray")
