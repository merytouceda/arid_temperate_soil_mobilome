# plasmid analysis for the base dataset

library(tidyverse)
library(vegan)
library(see)
library(lme4)
library(car)
library(performance)
library(lubridate)
library(see)

pal1 <- c("#D81B60", "#FFC107", "#004D40")
pal2 <- c("#D81B60", "#997404", "#004D40")

#plasmids_counts <- read.table("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/plasmids_count_table.txt", sep = "\t", header = T)
# with PTUs 
plasmids_counts <- read.table("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/ptus_count_table.txt", sep = "\t", header = T)
#md <- read.table("/Volumes/BunnyBike/mge_urban/base/base_metadata.txt", sep = "\t", header = T)
md <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_base_march25.csv")

plasmids_counts_long <- plasmids_counts %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "Counts") %>%
  mutate(Sample = gsub(Sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(Sample = gsub(Sample, pattern = "X", replacement = ""))

md_filtered <- md %>%
  mutate(Sample = gsub(sample_id, pattern = "_[0-9]", replacement = "")) %>%
  filter(Sample %in% plasmids_counts_long$Sample) %>%
  filter(!general_env_feature == "Australia") %>%
  filter(!is.na(vegetation_type)) %>%
  filter(!vegetation_type == "Forest")%>%
  filter(!Sample == "8459")

plasmids_counts_filtered <- plasmids_counts_long %>%
  filter(Sample %in% md_filtered$Sample) %>%
  pivot_wider(names_from = Contig, values_from = Counts ) %>%
  arrange(desc(Sample)) %>%
  column_to_rownames(var = "Sample") 

md_final<- md_filtered %>%
  distinct(Sample, .keep_all = T) %>%
  arrange(desc(Sample)) %>%   # check that the
  column_to_rownames(var = "Sample") 


# Read the final final metadata (the one with total reads)
md_final_final <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_march25.csv") %>%
  filter(!Sample == "8459") %>%
  column_to_rownames(var = "Sample") %>%
  # Add season
  mutate(
    # First convert the string dates to proper date objects
    date_clean = as_date(str_remove(collection_date, " UTC")),
    # Extract month number
    month = month(date_clean),
    # Create season based on meteorological seasons (Southern Hemisphere)
    season = case_when(
      month %in% c(12, 1, 2) ~ "Summer",
      month %in% c(3, 4, 5) ~ "Fall",
      month %in% c(6, 7, 8) ~ "Winter",
      month %in% c(9, 10, 11) ~ "Spring",
      TRUE ~ NA_character_
    )
  )


# make sure the rows are in the same order
rownames(plasmids_counts_filtered) == rownames(md_final)
rownames(plasmids_counts_filtered) == rownames(md_final_final) # final final

summary(rowSums(plasmids_counts_filtered))

# ------------------------------------------------------------------------------------------------------------------------------
# ALPHA DIVERSITY
###########

md_final_final$plasmid_rich <- specnumber(rrarefy(plasmids_counts_filtered, sample = 3101))
md_final_final$plasmid_shannon <- diversity(rrarefy(plasmids_counts_filtered, sample = 3101))

md_final_final$combined_var <- paste(md_final_final$general_env_feature, md_final_final$vegetation_type)

ggplot(md_final_final, aes(x = season, y = plasmid_rich, color= combined_var))+
  #geom_violin(aes(fill =urban.natural, color = urban.natural), binwidth = 1, alpha=0.4)+
  geom_boxplot(alpha=0.6, outlier.shape = NA)+
  #geom_point(aes(), size = 1, alpha = 0.7)+
  xlab(NULL) + 
  ylab("Plasmid richness")+
  scale_fill_see()+
  scale_color_see()+
  theme_bw()+
  theme(legend.position = "bottom", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
  #facet_wrap(~general_env_feature)
ggsave("/Volumes/BunnyBike/mge/base/figures/plasmid_rich.pdf", device = "pdf", width = 5, height = 5 , units = "in")


md_final_final$combined_var <- factor(md_final_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                                 "Arid Shrubland", "Temperate Shrubland", 
                                                                                 "Arid Woodland", "Temperate Woodland"))
# nice plot
ggplot(md_final_final, aes(x = vegetation_type, y = plasmid_shannon, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  ylab("Plasmid richness")+ 
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_shannon_paralell.pdf", device = "pdf", width = 5, height = 5 , units = "in")


Anova(lmer(plasmid_rich~  vegetation_type * general_env_feature + (1|season), data = md_final_final), type = 3)
#vegetation_type                        1.1912  2    0.55124    
#general_env_feature                    0.7315  1    0.39240    
#vegetation_type:general_env_feature    6.5910  2    0.03705 * 
r2(lmer(plasmid_rich ~  vegetation_type * general_env_feature + (1|season), data = md_final_final))


Anova(lmer(plasmid_shannon~  vegetation_type * general_env_feature + (1|season), data = md_final_final), type = 3)
#vegetation_type                        0.4867  2    0.78398    
#general_env_feature                    0.0002  1    0.98999    
#vegetation_type:general_env_feature    7.8109  2    0.02013 *
r2(lmer(plasmid_shannon ~  vegetation_type * general_env_feature + (1|season), data = md_final_final))

colnames(md_final_final)
# richness as number of plasmids/number of reads (Albert's idea)
'''
md_final <- md_final %>%
  mutate(plasmids_rich_2 = rowSums(plasmids_counts_filtered)/tot_reads)

ggplot(md_final, aes(x = vegetation_type, y = plasmids_rich_2, color = vegetation_type, fill = vegetation_type))+
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
  ylab("Plasmid diversity")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_rich_2_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")


anova(lm(plasmids_rich_2 ~  vegetation_type * general_env_feature, data = md_final))
'''

# ------------------------------------------------------------------------------------------------------------------------------
##############
# RPKM
##############
# 1. prepare the length and the total reads to left join
# import and format the table of lengths per gene
all_summary_filtered <- read_csv("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/circular_10kb_plasmids.csv")

plasmids_len <- all_summary_filtered %>%
  dplyr::select(c("seq_name", "length"))

colnames(plasmids_len) <- c("Contig", "len")

# select the total reads per sample from the metadata
total_reads <- md_final_final %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(Sample, tot_reads)

# 2. pivot longer the count table to make every gene in every sample a row
plasmids_long <- plasmids_counts_filtered%>%
  rownames_to_column(var = "Sample") %>%
  pivot_longer(cols = -c("Sample"), values_to = "Count", names_to= "Contig") 

# 3. join the total reads to each row based on sample
plasmids_long_plusreads <- left_join(plasmids_long, total_reads, by = "Sample")
# 4. join the gene length to each row based on gene
plasmids_long_plusreads_plustlen <- left_join(plasmids_long_plusreads, plasmids_len, by = "Contig")

# 5. compute the RPKM calculation and create count table from it 
rpkm_plasmids_final <- plasmids_long_plusreads_plustlen  %>%
  dplyr::mutate(RPKM = (Count*10e6)/(tot_reads*as.numeric(len))) %>%
  dplyr::select(-c("Count", "tot_reads", "len")) %>%
 # arrange(Sample)%>% # use this for the mantel tests
  pivot_wider(names_from = Sample, values_from = RPKM)%>%
  column_to_rownames(var = "Contig") %>%
  unique() %>%
  drop_na()


# ------------------------------------------------------------------------------------------------------------------------------
############################
# beta-diversity
############################
plasmids.bray <- vegdist(t(rpkm_plasmids_final), method="bray")

#ordination (non-multidimensional scaling)
plasmids.nmds <- metaMDS(plasmids.bray, k=2, try = 100)
md_final_final$Axis01 = plasmids.nmds$points[,1]
md_final_final$Axis02 = plasmids.nmds$points[,2]
plasmids.nmds$stress #0.03989818 the smaller the better, good <0.3)


md_final_final$combined_var <- paste(md_final_final$general_env_feature, md_final_final$vegetation_type)
md_final_final$combined_var <- factor(md_final_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                  "Arid Shrubland", "Temperate Shrubland", 
                                                                  "Arid Woodland", "Temperate Woodland"))
ggplot(md_final_final, aes(Axis01, Axis02))+
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
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_beta_rpkm_combinedvar.pdf", device = "pdf", width = 5, height = 5 , units = "in")

adonis2(plasmids.bray ~ season + vegetation_type * general_env_feature, data = md_final_final,  permutations = 999, method = "bray")

adonis2(plasmids.bray ~ vegetation_type * general_env_feature, data = md_final_final,  
        permutations = 999, method = "bray", strata = md_final_final$month)

# check other variables
adonis2(plasmids.bray ~ vegetation_type + general_env_feature + organic_carbon + exc_magnesium + ph,
        data = md_final_final,  permutations = 999, method = "bray")
#organic_carbon        1   1.2177 0.04181  8.6905  0.001 ***
#  exc_magnesium         1   0.9745 0.03346  6.9550  0.001 ***
#  ph                    1   0.9963 0.03421  7.1105  0.001 ***




# for mantel tests
md_arid <- md_final %>%
  filter(general_env_feature == "Arid") %>%
  mutate(sample = gsub(sample_id, pattern = "_[0-9]$", replacement = ""))

rpkm_plasmids_final_arid <- rpkm_plasmids_final %>% 
  select(md_arid$sample)

rpkm_plasmids_final_temperate <- rpkm_plasmids_final %>% 
  select(!md_arid$sample)

plasmids.bray.arid <- vegdist(t(rpkm_plasmids_final_arid), method="bray")
plasmids.bray.temperate <- vegdist(t(rpkm_plasmids_final_temperate), method="bray")



# ------------------------------------------------------------------------------------------------------------------------------
############################
# abundance
############################
sum_rel_ab_plasmids <- as.data.frame(t(rpkm_plasmids_final)) 

sum_rel_ab_plasmids_final <- sum_rel_ab_plasmids %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(plasmids_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(plasmids_relab) %>%
  rownames_to_column(var = "Sample")

md_final_final<- md_final_final %>%
  rownames_to_column(var = "Sample") %>%
  left_join(., sum_rel_ab_plasmids_final, by = "Sample")


ggplot(md_final_final, aes(x = vegetation_type, y = plasmids_relab))+
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
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/plasmid_relab_veg.pdf", device = "pdf", width = 5, height = 5 , units = "in")


ggplot(md_final_final, aes(x = vegetation_type, y = plasmids_relab, color = vegetation_type, fill = vegetation_type))+
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
  ylab("Plasmid relative abundance")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_relab_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")


anova(lm(plasmids_relab ~  vegetation_type * general_env_feature, data = md_final))

Anova(lmer(plasmids_relab ~  vegetation_type * general_env_feature  + (1|season), data = md_final_final), type = 3)
#vegetation_type                      1.1669  2    0.55797    
#general_env_feature                  0.3043  1    0.58123    
#vegetation_type:general_env_feature  6.3476  2    0.04184 * 



# ------------------------------------------------------------------------------------------------------------------------------
############################
# length average plot
############################

counts_pab <- rpkm_plasmids_final
counts_pab[counts_pab > 0] <- 1
md_final_final_join <- md_final_final #%>%
  rownames_to_column(var = "Sample")

counts_pab_length <- counts_pab %>%
  ungroup() %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(!Contig, names_to = "Sample", values_to = "pab") %>%
  left_join(plasmids_len, by = "Contig") %>%
  filter(pab == 1) %>%
  dplyr::group_by(Sample) %>%
  dplyr::summarize(av_len = mean(len, na.rm = TRUE), 
            min_len = min(len, na.rm = TRUE), 
            max_len = max(len, na.rm = TRUE)) %>%
  left_join(md_final_final_join, by = "Sample")


ggplot(counts_pab_length, aes(x = vegetation_type, y = av_len))+
  #geom_violin(aes(fill =urban.natural, color = urban.natural), binwidth = 1, alpha=0.4)+
  geom_boxplot(alpha=0.6, outlier.shape = NA )+
  geom_jitter(aes(color= vegetation_type), size = 2.5)+
  xlab(NULL) + 
  ylab("Plasmid average length")+
  scale_fill_see()+
  scale_color_see()+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~general_env_feature)
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/plasmid_avlength_veg.pdf", device = "pdf", width = 5, height = 5 , units = "in")


ggplot(counts_pab_length, aes(x = vegetation_type, y = av_len, color = vegetation_type, fill = vegetation_type))+
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
  ylab("Plasmid length")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_length_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")


anova(lm(av_len ~  vegetation_type * general_env_feature, data = counts_pab_length))

Anova(lmer(av_len ~  vegetation_type * general_env_feature  + (1|season), data = counts_pab_length), type = 3)
#vegetation_type                         3.5443  2    0.16997    
#general_env_feature                     0.7058  1    0.40084    
#vegetation_type:general_env_feature     7.3741  2    0.02505 * 
r2(lmer(av_len ~  vegetation_type * general_env_feature  + (1|season), data = counts_pab_length))


# new aes
# nice plot
counts_pab_length$combined_var <- paste(counts_pab_length$general_env_feature, counts_pab_length$vegetation_type)
counts_pab_length$combined_var <- factor(counts_pab_length$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                              "Arid Shrubland", "Temperate Shrubland", 
                                                              "Arid Woodland", "Temperate Woodland"))
ggplot(counts_pab_length, aes(x = vegetation_type, y = av_len, color = combined_var, fill = combined_var))+
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
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_length_nice_parallel.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# ------------------------------------------------------------------------------------------------------------------------------
# length frequency distributions

counts_pab <- rpkm_plasmids_final
counts_pab[counts_pab > 0] <- 1

md_final_final_join <- md_final_final %>%
  rownames_to_column(var = "Sample") %>%
  dplyr::select(c("Sample", "general_env_feature", "vegetation_type"))

counts_pab_length <- counts_pab %>%
  ungroup() %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(!Contig, names_to = "Sample", values_to = "pab") %>%
  left_join(plasmids_len, by = "Contig") %>%
  filter(pab == 1)%>%
  left_join(md_final_final_join, by = "Sample")



ggplot(counts_pab_length, aes(x = len, fill = vegetation_type))+
  geom_histogram(position = "identity", alpha = 0.6, bins = 50) +
  facet_wrap(~general_env_feature)+
  scale_fill_manual(values = pal2)+
  ylab("Frequency") + 
  xlab("Plasmid length")+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_length_distribution.pdf", device = "pdf", width = 7, height = 5 , units = "in")





# ------------------------------------------------------------------------------------------------------------------------------
## Tukey test & visualizing the interactions
library(emmeans)
library(multcomp)

m1 <- lmer(plasmid_shannon~  vegetation_type * general_env_feature  + (1|season), data = md_final_final)
m2 <- lmer(plasmids_relab ~  vegetation_type * general_env_feature  + (1|season), data = md_final_final)
m3 <- lmer(av_len ~  vegetation_type * general_env_feature  + (1|season), data = counts_pab_length)


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
       y = "Plasmid Length") +
  theme_bw() +
  theme(legend.position = "none",
        #axis.text.x = element_text(angle = 45, hjust = 1),
        panel.grid.minor = element_blank(), text = element_text(size = 16)) +
  scale_color_manual(values = pal1)
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ptus_length_interaction.pdf", device = "pdf", width = 3, height = 5 , units = "in")



# ------------------------------------------------------------------------------------------------------------------------------
# save md_final for correlations
# add the average plasmid length
av_len <- counts_pab_length %>%
  select(c("Sample", "av_len"))

md_final_final <- md_final_final %>%
  left_join(.,av_len, by = "Sample")

write_csv(md_final_final, "/Volumes/BunnyBike/mge_urban/base/md_final_final_ptus.csv")



############################
# beta-diversity separated arid and temperate
############################
md_arid <- md_final_final %>%
  rownames_to_column(var = "Sample") %>%
  filter(general_env_feature == "Arid")

rpkm_plasmids_arid <- rpkm_plasmids_final %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "Count") %>%
  filter(Sample %in% md_arid$Sample) %>%
  pivot_wider(names_from = "Sample", values_from = "Count") %>%
  column_to_rownames(var = "Contig")

rpkm_plasmids_temperate <- rpkm_plasmids_final %>%
  rownames_to_column(var = "Contig") %>%
  pivot_longer(-Contig, names_to = "Sample", values_to = "Count") %>%
  filter(!Sample %in% md_arid$Sample) %>%
  pivot_wider(names_from = "Sample", values_from = "Count") %>%
  column_to_rownames(var = "Contig")

plasmids.bray_arid <- vegdist(t(rpkm_plasmids_arid), method="bray")
plasmids.bray_temperate <- vegdist(t(rpkm_plasmids_temperate), method="bray")


