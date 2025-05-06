# analysis of the bacterial traits in base project

library(tidyverse)
library(vegan)
library(ggpubr)
library(performance)
library(car)
library(lme4)
pal1 <- c("#D81B60", "#FFC107", "#004D40")


#md <- read.table("/Volumes/BunnyBike/mge_urban/base/base_metadata.txt", sep = "\t", header = T)
md <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_bact.csv")

md_final <- md %>%
  mutate(Sample = gsub(sample_id, pattern = "_[0-9]", replacement = "")) %>%
  filter(!general_env_feature == "Australia") %>%
  filter(!is.na(vegetation_type)) %>%
  filter(!vegetation_type == "Forest") %>%
  arrange(desc(Sample)) %>%
  filter(!Sample == "8459")

# ------------------------------------------------------------------------------------------------------------------------------
# Average genome size (AGS)
ags_result <- read.table("/Volumes/BunnyBike/mge_urban/base/traits/base_ags_result.txt", sep = "\t", header = T)
colnames(ags_result) <- c("Sample", "totalbp", "genome_number", "ags")

ags_result <- ags_result %>%
  mutate(Sample = as.character(Sample))

md_ags <- md_final %>%
  filter(Sample %in% ags_result$Sample) %>%
  left_join(., ags_result, by = "Sample") 

## I am missing samples in the metadata, why? 189 samples total but when I join with metadata I have 143 (I think those are the australia, forest, etc ones)

ggplot(md_ags, aes(x = vegetation_type, y = ags))+
  geom_boxplot(aes(color = vegetation_type))+
  geom_jitter(aes(color = vegetation_type))+
  facet_wrap(~general_env_feature)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average genome size (AGS)")+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ags.pdf", device = "pdf", width = 8, height = 8 , units = "in")

anova(lm(ags ~  vegetation_type * general_env_feature, data = md_ags))

# nice plot
md_ags$combined_var <- paste(md_ags$general_env_feature, md_ags$vegetation_type)
md_ags$combined_var <- factor(md_ags$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                  "Arid Shrubland", "Temperate Shrubland", 
                                                                  "Arid Woodland", "Temperate Woodland"))
ggplot(md_ags, aes(x = vegetation_type, y = ags, color = combined_var, fill = combined_var))+
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
  ylab("Average genome size")+
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/ags_nice_parallel.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(ags ~  vegetation_type * general_env_feature + (1|season), data = md_ags), type = 3)
#vegetation_type                        0.5478  2     0.7604    
#general_env_feature                   16.4762  1  4.926e-05 ***
#vegetation_type:general_env_feature    2.4993  2     0.2866
r2(lmer(ags ~  vegetation_type * general_env_feature + (1|season), data = md_ags))

# ------------------------------------------------------------------------------------------------------------------------------
# Average copy number (ACN)
coverage_result <- read.table("/Volumes/BunnyBike/mge_urban/base/traits/base__16S_coverage.txt", sep = "\t", header = T)
colnames(coverage_result) <- c("Sample", "coverage")

coverage_result <- coverage_result %>%
  mutate(Sample = as.character(Sample))

md_acn <- md_ags %>%
  filter(Sample %in% coverage_result$Sample) %>%
  left_join(., coverage_result, by = "Sample") %>%
  mutate(acn = coverage/genome_number)

## I am missing samples in the metadata, why? 189 samples total but when I join with metadata I have 143 (I think those are the australia, forest, etc ones)

ggplot(md_acn, aes(x = vegetation_type, y = acn))+
  geom_boxplot(aes(color = vegetation_type))+
  geom_jitter(aes(color = vegetation_type))+
  facet_wrap(~general_env_feature)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average copy number (ACN)")+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/acn.pdf", device = "pdf", width = 8, height = 8 , units = "in")

anova(lm(acn ~  vegetation_type * general_env_feature, data = md_acn))

# nice plot
ggplot(md_acn, aes(x = vegetation_type, y = acn, color = combined_var, fill = combined_var))+
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
  ylab("Average copy number (ACN)")+
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/acn_nice_parallel.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(acn ~  vegetation_type * general_env_feature + (1|season), data = md_acn), type = 3)
#vegetation_type                        9.8831  2  0.0071436 ** 
#general_env_feature                   13.5186  1  0.0002362 ***
#vegetation_type:general_env_feature    3.8702  2  0.1444065 
r2(lmer(acn ~  vegetation_type * general_env_feature + (1|season), data = md_acn))

# ------------------------------------------------------------------------------------------------------------------------------
# GC content
gc_result <- read.table("/Volumes/BunnyBike/mge_urban/base/traits/base_gc_result.txt", sep = "\t", header = T)

gc_result  <- gc_result  %>%
  mutate(Sample = as.character(Sample))

md_gc <- md_acn %>%
  filter(Sample %in% gc_result$Sample) %>%
  left_join(., gc_result, by = "Sample") 


ggplot(md_gc, aes(x = vegetation_type, y = gc_mean))+
  geom_boxplot(aes(color = vegetation_type))+
  geom_jitter(aes(color = vegetation_type))+
  facet_wrap(~general_env_feature)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average GC content (%)")+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/gc_mean.pdf", device = "pdf", width = 8, height = 8 , units = "in")

anova(lm(gc_mean ~  vegetation_type * general_env_feature, data = md_gc))


ggplot(md_gc, aes(x = vegetation_type, y = gc_var))+
  geom_boxplot(aes(color = vegetation_type))+
  geom_jitter(aes(color = vegetation_type))+
  facet_wrap(~general_env_feature)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Average GC variance (%)")+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/gc_var.pdf", device = "pdf", width = 8, height = 8 , units = "in")

anova(lm(gc_var~  vegetation_type * general_env_feature, data = md_gc))


# nice plot
ggplot(md_gc, aes(x = vegetation_type, y = gc_mean, color = vegetation_type, fill = vegetation_type))+
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
  ylab("GC content (%)")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/gc_mean_nice.pdf", device = "pdf", width = 7, height = 5 , units = "in")

Anova(lmer(gc_mean ~  vegetation_type * general_env_feature + (1|season), data = md_gc), type = 3)
#vegetation_type                         1.7109  2     0.4251    
#general_env_feature                    49.3067  1  2.189e-12 ***
#vegetation_type:general_env_feature     3.4661  2     0.1767  
r2(lmer(gc_mean ~  vegetation_type * general_env_feature + (1|season), data = md_gc))

ggplot(md_gc, aes(x = vegetation_type, y = gc_mean, color = combined_var, fill = combined_var))+
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
  ylab("GC content (%)")+
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/gc_mean_nice_parallel.pdf", device = "pdf", width = 5, height = 5 , units = "in")




write_csv(md_gc, "/Volumes/BunnyBike/mge_urban/base/md_final_final_traits.csv")
# ------------------------------------------------------------------------------------------------------------------------------
# correlations between traits


# AGS vs ACN
ggplot(md_gc, aes(x = ags, y = acn))+
  geom_point()+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Average genome size (AGS)")+
  ylab("Average copy number (ACN)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ags_acn_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# AGS vs GC content
ggplot(md_gc, aes(x = ags, y = gc_mean))+
  geom_point()+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Average genome size (AGS)")+
  ylab("Average GC content (%)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ags_gc_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# ACN vs GC content
ggplot(md_gc, aes(x = acn, y = gc_mean))+
  geom_point()+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Average copy number (ACN)")+
  ylab("Average GC content (%)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_acn_gc_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# ------------------------------------------------------------------------------------------------------------------------------
## Tukey test & visualizing the interactions
library(emmeans)
library(multcomp)

m1 <- lmer(ags~  vegetation_type * general_env_feature  + (1|season), data = md_gc)
m2 <- lmer(acn ~  vegetation_type * general_env_feature  + (1|season), data = md_gc)
m3 <- lmer(gc_mean~  vegetation_type * general_env_feature  + (1|season), data = md_gc)

# Step 1: Get estimated marginal means for the interaction
emm_interaction <- emmeans(m3, ~ vegetation_type * general_env_feature)

# Step 2: Perform pairwise comparisons with Tukey adjustment
pairs_result <- pairs(emm_interaction, adjust = "tukey")

# Step 3: View the results
summary(pairs_result)


# Optional: Create a compact letter display
cld_result <- cld(emm_interaction, Letters = letters, adjust = "tukey")
cld_result





