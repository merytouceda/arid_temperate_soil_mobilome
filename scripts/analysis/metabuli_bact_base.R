# parse metabuli output base bacteria

library(tidyverse)
library(car)
library(lme4)
library(car)
library(performance)
library(vegan)
pal1 <- c("#D81B60", "#FFC107", "#004D40")

# load wd
setwd("/Volumes/BunnyBike/mge_urban/base/metabuli/contigs/reports/parsed/")


# initialize empty table
metabuli_count <- data.frame(matrix(NA, nrow = 0, ncol = 10))  #after reviewing this I think the "wide" format is most useful since it keep the LOD information
colnames(metabuli_count) <- c("reads", "taxon", "superkingdom" ,"phylum" ,"class" ,"order","family" , "genus" , "species", "sample")

# Retrieve all samples info and add to this empty table
files <- list.files(pattern="*.tsv", full.names=F, recursive=FALSE)

for (file in files){
  # save file name
  filename <- basename(file)
  
  #open file (raise error if it doesn't open)
  d <- tryCatch(read_tsv(file), error=function(e) NULL)
  
  taxonomy_table <- d %>%
    dplyr::select(c("reads","taxon", "superkingdom" ,"phylum" ,"class" ,"order","family" , "genus" , "species")) %>%
    mutate(sample = filename) %>%
    mutate_if(is.numeric, ~replace_na(.x, 0))
  
  # append to form a big dataset
  metabuli_count <-rbind(metabuli_count, taxonomy_table)
  
}


# create count table (samples as columns, taxa as rows)
metabuli_count_table <- metabuli_count %>%
  dplyr::mutate(sample = gsub(sample, pattern = "base_contigs_", replacement = "")) %>%
  dplyr::mutate(sample = gsub(sample, pattern = ".parsed.tsv", replacement = "")) %>%
  pivot_wider(names_from = "sample", values_from = "reads") %>%
  mutate_if(is.numeric, ~replace_na(.x, 0))


# calculate percentage unannotated
unk_count_table <- metabuli_count_table %>%
  filter(taxon %in% c("unclassified", "root", "Bacteria"))

# calculate percentages
mean(colSums(unk_count_table[,9:197]) / colSums(metabuli_count_table[,9:197]))
#  0.7128641
sd(colSums(unk_count_table[,9:197]) / colSums(metabuli_count_table[,9:197]))
# 0.01028357


# remove anything that is not a classified bacteria or archaea
final_count_table <- metabuli_count_table %>%
  filter(!is.na(phylum)) %>%
  filter(!superkingdom %in% c("Viruses", "Eukaryota"))

# retrieve viruses just for funsies
virus_count_table <- metabuli_count_table %>%
  filter(superkingdom == "Viruses")

# retrieve euk
euk_count_table <- metabuli_count_table %>%
  filter(superkingdom == "Eukaryota")


# calculate richness

# find duplicates
length(unique(final_count_table$taxon))
duple <- final_count_table[duplicated(final_count_table$taxon),]

# remove duplicates ("environmental samples") by assigning them the genus and summarizing 
taxon_count_table <- final_count_table %>%
  mutate(taxon = case_when(taxon == "environmental samples" ~ genus, 
                           TRUE ~taxon)) %>%
  dplyr::select(-c("superkingdom", "phylum", "class", "order", "family", "genus", "species")) %>%
  group_by(taxon) %>%
  summarize_all(sum)%>%
  column_to_rownames(var = "taxon")

# filter taxon table
# Prevalence
# Keep those columns (species) that are prevalent
# meaning, their abundance has to be > 0 in at least 10% of the samples (19 samples)
# get the taxons that have the required prevalence
taxon_count_table_prvalence <- taxon_count_table %>%
  rownames_to_column(var = "taxon") %>%
  pivot_longer(-taxon, names_to = "sample", values_to = "count") %>%
  mutate(pab = case_when(count == 0 ~ 0, 
                         TRUE ~ 1)) %>%
  group_by(taxon) %>%
  summarize(pab_total = sum(pab)) %>%
  filter(pab_total > 19)

#filter the table
filter_taxon_count_table <- taxon_count_table %>%
  rownames_to_column(var = "taxon") %>%
  filter(taxon %in% taxon_count_table_prvalence$taxon)

# Abundance  ?? ask ALISE
# calculate summary of 
filter_taxon_count_table_long <- filter_taxon_count_table %>%
  pivot_longer(-taxon, names_to = "sample", values_to = "counts") %>%
  group_by(taxon) %>%
  summarize(suma = sum(counts)) %>%
  filter(suma > 100)

filter_taxon_count_table_final <- filter_taxon_count_table %>%
  filter(taxon %in% filter_taxon_count_table_long$taxon)

summary(filter_taxon_count_table_long$`sum(counts)`)
hist(filter_taxon_count_table_long$`sum(counts)`)

#write_csv(filter_taxon_count_table_final, "/Volumes/BunnyBike/mge_urban/base/metabuli/metabuli_count_table_final.csv")
# ----------------------------------------------------------------------------------------------------------------
##### RICHNESS
# load metadata
md <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_votus.csv")

filter_taxon_count_table_long <- filter_taxon_count_table_final %>%
  pivot_longer(-taxon, names_to = "sample", values_to = "count")

# filter md
md_filtered <- md %>%
  #mutate(sample = gsub(sample_id, pattern = "_[0-9]", replacement = "")) %>%
  mutate(sample = Sample) %>%
  filter(sample %in% filter_taxon_count_table_long$sample) %>%
  filter(!general_env_feature == "Australia") %>%
  filter(!is.na(vegetation_type)) %>%
  filter(!vegetation_type == "Forest") %>%
  filter(!sample == "8459") # the sample with no collection date and thus no season

filter_taxon_count_table_fitered <- filter_taxon_count_table_long %>%
  filter(sample %in% md_filtered$sample) %>%
  pivot_wider(names_from = taxon, values_from = count) %>%
  arrange(desc(as.numeric(sample))) %>%
  column_to_rownames(var = "sample") 

md_final<- md_filtered %>%
  arrange(desc(sample)) %>%   # check that the
  column_to_rownames(var = "sample") 


# make sure the rows are in the same order
rownames(filter_taxon_count_table_fitered) == rownames(md_final)


# Check what number to rarify by
summary(rowSums(filter_taxon_count_table_fitered))

# Richness
md_final$bact_rich <- specnumber(rrarefy(filter_taxon_count_table_fitered, sample = 62443))
# SHANNON with rarefaction
md_final$bact_shannon <- diversity(rrarefy(filter_taxon_count_table_fitered, sample = 62443), index = "shannon")


# nice plot
ggplot(md_final, aes(x = vegetation_type, y = bact_rich, color = vegetation_type, fill = vegetation_type))+
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
  ylab("Bacterial richness")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/bact_rich_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")

Anova(lmer(bact_rich ~  vegetation_type * general_env_feature + (1|season), data = md_final), type = 3)
Anova(lmer(bact_rich ~   general_env_feature * vegetation_type  + (1|season), data = md_final), type = 3)
#vegetation_type                        4.3868  2     0.1115    
#general_env_feature                   26.6924  1  2.386e-07 ***
#vegetation_type:general_env_feature    2.4284  2     0.2970  

library(performance)
r2(lmer(bact_rich ~  vegetation_type * general_env_feature + (1|season), data = md_final))

Anova(lmer(bact_shannon ~  vegetation_type * general_env_feature + (1|season), data = md_final), type = 3)
#vegetation_type                        3.9640  2     0.1378    
#general_env_feature                    4.1516  1     0.0416 *  
#vegetation_type:general_env_feature    1.3855  2     0.5002 
# alternative aesthetic
r2(lmer(bact_shannon~  vegetation_type * general_env_feature + (1|season), data = md_final))


md_final$combined_var <- paste(md_final$general_env_feature, md_final$vegetation_type)
md_final$combined_var <- factor(md_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                  "Arid Shrubland", "Temperate Shrubland", 
                                                                  "Arid Woodland", "Temperate Woodland"))

ggplot(md_final, aes(x = vegetation_type, y = bact_rich, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  ylab("Bacterial Shannon diversity")+ 
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/bact_shannon_paralell.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# ----------------------------------------------------------------------------------------------------------------
##############
# RPM
##############
# select the total reads per sample from the metadata
total_reads <- md_final %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(sample, tot_reads)

# 2. pivot longer the count table to make every gene in every sample a row
filter_taxon_count_table_fitered_long <- filter_taxon_count_table_fitered %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(cols = -c("sample"), values_to = "Count", names_to= "taxon") 

# 3. join the total reads to each row based on sample
filter_taxon_count_table_fitered_long_plusreads <- left_join(filter_taxon_count_table_fitered_long, total_reads, by = "sample")


# 5. compute the RPKM calculation and create count table from it 
rpm_bact_final <- filter_taxon_count_table_fitered_long_plusreads  %>%
  dplyr::mutate(RPM = (Count*10e6)/(tot_reads)) %>% # remove the 106?
  dplyr::select(-c("Count", "tot_reads")) %>%
  #arrange(sample)%>% # use only for Mantel tests
  pivot_wider(names_from = sample, values_from = RPM)%>%
  column_to_rownames(var = "taxon") %>%
  unique() %>%
  drop_na()


# ----------------------------------------------------------------------------------------------------------------
############################
# beta-diversity
############################

bact.bray <- vegdist(t(rpm_bact_final), method="bray")

#ordination (non-multidimensional scaling)
bact.nmds <- metaMDS(bact.bray , k=2, try = 100)
md_final$Axis01 = bact.nmds$points[,1]
md_final$Axis02 = bact.nmds$points[,2]
bact.nmds$stress #0.03989818 the smaller the better, good <0.3)

colnames(md_final)
ggplot(md_final, aes(Axis01, Axis02))+
  geom_point(aes(color= vegetation_type, shape = general_env_feature), size=4)+
  #stat_ellipse(aes(color = vegetation_type))+
  stat_ellipse(aes(group = general_env_feature, linetype = general_env_feature))+
  scale_shape_manual(values=c(16, 13))+
  scale_color_manual(values = pal1)+
  theme_bw()+
  theme(legend.position="none", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/bact_communitycomp.pdf", device = "pdf", width = 7, height = 5 , units = "in")

adonis2(bact.bray ~ season + vegetation_type * general_env_feature, data = md_final,  permutations = 999, method = "bray")
#this is the way

adonis2(bact.bray ~ vegetation_type + general_env_feature + ph  + organic_carbon , data = md_final,  permutations = 999, method = "bray")
#ph                    1   1.2760 0.06278 16.929  0.001 ***
#organic_carbon        1   1.0777 0.05303 14.299  0.001 ***

md_final$combined_var <- paste(md_final$general_env_feature, md_final$vegetation_type)
md_final$combined_var <- factor(md_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                  "Arid Shrubland", "Temperate Shrubland", 
                                                                  "Arid Woodland", "Temperate Woodland"))

# new aesthetic
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
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/bact_beta_rpm_combinedvar.pdf", device = "pdf", width = 5, height = 5 , units = "in")




# for mantel tests
md_arid <- md_final %>%
  filter(general_env_feature == "Arid") %>%
  mutate(sample = gsub(sample_id, pattern = "_[0-9]$", replacement = ""))

rpm_bact_final_arid <- rpm_bact_final %>% 
  select(md_arid$sample)

rpm_bact_final_temperate <- rpm_bact_final %>% 
  select(!md_arid$sample)

bact.bray.arid <- vegdist(t(rpm_bact_final_arid), method="bray")
bact.bray.temperate <- vegdist(t(rpm_bact_final_temperate), method="bray")


# ----------------------------------------------------------------------------------------------------------------
############################
# abundance
############################
sum_rel_ab_bact <- as.data.frame(t(rpm_bact_final)) 

sum_rel_ab_bact_final <- sum_rel_ab_bact %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(bact_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(bact_relab) %>%
  rownames_to_column(var = "Sample")

md_final<- md_final %>%
  rownames_to_column(var = "Sample") %>%
  left_join(., sum_rel_ab_bact_final, by = "Sample")


ggplot(md_final, aes(x = vegetation_type, y = bact_relab, color = vegetation_type, fill = vegetation_type))+
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
  ylab("Bacterial relative abundance")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/bact_relab_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")


anova(lm(bact_relab ~  vegetation_type * general_env_feature, data = md_final))
#Response: bact_relab
#Df     Sum Sq    Mean Sq F value  Pr(>F)  
#vegetation_type                       2 6.3645e+09 3182274493  1.5440 0.21722  
#general_env_feature                   1 3.6696e+09 3669602612  1.7804 0.18431  
#vegetation_type:general_env_feature   2 1.1044e+10 5521825654  2.6791 0.07222 .
#Residuals                           137 2.8237e+11 2061089654  


write_csv(md_final, "/Volumes/BunnyBike/mge_urban/base/md_final_bact.csv")




# ------------------------------------------------------------------------------------------------------------------------------
## Tukey test & visualizing the interactions
library(emmeans)
library(multcomp)

m1 <- lmer(bact_rich~  general_env_feature *vegetation_type  + (1|season), data = md_final)
m2 <- lmer(bact_shannon~  vegetation_type * general_env_feature  + (1|season), data = md_final)

#library(modelbased)
#estimate_contrasts(m2, contrast = vegetation_type)

# Step 1: Get estimated marginal means for the interaction
emm_interaction <- emmeans(m2, ~ vegetation_type * general_env_feature)

# Step 2: Perform pairwise comparisons with Tukey adjustment
pairs_result <- pairs(emm_interaction, adjust = "tukey")

# Step 3: View the results
summary(pairs_result)


# Optional: Create a compact letter display
cld_result <- cld(emm_interaction, Letters = letters, adjust = "tukey")
cld_result









# ----------------------------------------------------------------------------------------------------------------
# taxonomy


# make rpm into presence_abscence table
#counts_pab <- rpm_bact_final
#counts_pab[counts_pab > 0] <- 1

md_join <- md_final %>%
  mutate(sample = as.character(Sample)) %>%
  dplyr::select(c("sample", "general_env_feature", "vegetation_type"))

# Join the pab table to the taxonomy
taxonomy_pab <- rpm_bact_final %>%
  rownames_to_column(var = "taxon") %>%
  left_join(., taxonomy_table, by = "taxon") %>%
  dplyr::select(-c("reads", "sample")) %>%
  dplyr::select(-c("taxon", "superkingdom", "class", "order", "family", "genus", "species")) %>%
  pivot_longer(-phylum, names_to = "sample", values_to = "pab") %>%
  left_join(., md_join, by = "sample") %>%
  group_by(phylum, general_env_feature, vegetation_type) %>%
  summarize(counts = sum(pab)) %>%
  mutate(pab = case_when(counts == 0 ~ 0, 
                         TRUE ~ 1)) %>%
  mutate(combined_var = paste(general_env_feature, vegetation_type)) %>%
  arrange(desc(counts)) %>%
  # calculate percentage to plot
  group_by(combined_var) %>%
  mutate(percentage = (counts / sum(counts)) * 100) %>%
  ungroup()

df_top10 <- taxonomy_pab %>%
  arrange(desc(percentage)) %>%
  mutate(rank = row_number()) %>%
  mutate(phylum = ifelse(rank > 49 | is.na(phylum), "Other", phylum)) %>%
  group_by(combined_var, phylum, general_env_feature, vegetation_type) %>%
  summarise(counts = sum(counts), percentage = sum(percentage), .groups = "drop") %>%
  arrange(combined_var, desc(percentage))

pal3 <- c("#FDAE61",  "#FEE08B", "#D53E4F" ,"#ABDDA4", 
          "purple","#9E0142","black" , "#F46D43","grey", "pink","#3288BD", "darkgreen")

ggplot(df_top10, aes(x = vegetation_type, y = percentage, fill = phylum, group = counts))+
  geom_col()+
  scale_fill_manual(values = pal3)+
  xlab("")+
  ylab("Proportion of counts (%)")+
  facet_wrap(~general_env_feature)+
  theme_bw()+
  theme(legend.position = "none")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/bact_phyla_nolegend.pdf", device = "pdf", width = 6, height = 5 , units = "in")





############################
# beta-diversity separated arid and temperate
############################
md_arid <- md_final_final %>%
  rownames_to_column(var = "Sample") %>%
  filter(general_env_feature == "Arid")

rpm_bact_arid <- rpm_bact_final %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "Count") %>%
  filter(Sample %in% md_arid$Sample) %>%
  pivot_wider(names_from = "Sample", values_from = "Count") %>%
  column_to_rownames(var = "Contig")

rpm_bact_temperate <- rpm_bact_final %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "Count") %>%
  filter(!Sample %in% md_arid$Sample) %>%
  pivot_wider(names_from = "Sample", values_from = "Count") %>%
  column_to_rownames(var = "Contig")

bact.bray_arid <- vegdist(t(rpm_bact_arid), method="bray")
bact.bray_temperate <- vegdist(t(rpm_bact_temperate ), method="bray")
