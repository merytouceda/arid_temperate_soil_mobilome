# correlation analyses between bacteria/plasmids/viruses for base project

library(tidyverse)
library(ggpubr)
pal1 <- c("#D81B60", "#FFC107", "#004D40")
pal1 <- c("#D81B60", "#997404", "#004D40")


#md_p <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_ptus.csv")
#md_v <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_votus.csv")
#md_b <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_bact.csv") # I have saved so that this already has all the variables

md_trait <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_traits.csv")
#md_virus_new <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_votus.csv") %>%
#  select(c("Sample", "ratio_vt"))
md_crispr <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_crispr.csv") %>%
  dplyr::select(c("Sample", "CRISPR_i"))
md_hostrange <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_hostrange.csv") %>%
  dplyr::select(c("sample", "prop_broad")) %>%
  dplyr::mutate(Sample = sample) %>%
  dplyr::select(-sample)

md_all <- md_trait %>%
  #left_join(md_virus_new, by = "Sample") %>%
  left_join(md_crispr, by = "Sample") %>%
  left_join(md_hostrange, by = "Sample")


###### INDIVIDUAL CORRELATIONS

# bacteria vs PTUS
ggplot(md_all, aes(x = bact_rich, y = plasmid_rich))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("Bacterial species richness")+
  ylab("PTU richness")+
  facet_wrap(~general_env_feature, scales = "free")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_pturich_bactrich.pdf", device = "pdf", width = 7, height = 5 , units = "in")

# final corr plot
ggplot(md_all, aes(x = bact_rich, y = plasmid_rich))+
  geom_point(aes(color = general_env_feature))+
  geom_smooth(method = "lm", aes(color = general_env_feature), se=F)+
  stat_cor(method = "pearson", aes(color = general_env_feature))+
  theme_bw() +
  xlab("Bacterial species richness")+
  ylab("PTU richness")+
  scale_color_manual(values = c("grey70", "grey30"))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_pturich_bactrich_final.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# suplementary figure making
# final corr plot
ggplot(md_all, aes(x = ags, y = gc_mean))+
  geom_point(aes(color = general_env_feature))+
  geom_smooth(method = "lm", aes(color = general_env_feature), se=F)+
  stat_cor(method = "pearson", aes(color = general_env_feature))+
  theme_bw() +
  xlab("AGS")+
  ylab("GC content")+
  scale_color_manual(values = c("grey70", "grey30"))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ags_gc_final.pdf", device = "pdf", width = 5, height = 5 , units = "in")









cor.test(md_all$bact_rich, md_all$plasmid_rich)

# by vegetation type
ggplot(md_all, aes(x = bact_rich, y = plasmid_rich))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type, fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Bacterial richness")+
  ylab("PTU richness")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_pturich_bactrich_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# bacteria vs vOTUS
ggplot(md_all, aes(x = bact_rich, y = virus_rich))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("Bacterial species richness")+
  ylab("vOTUS richness")+
  facet_wrap(~general_env_feature, scales = "free")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_voturich_bactrich.pdf", device = "pdf", width = 7, height = 5 , units = "in")

# final corr plot
ggplot(md_all, aes(x = bact_rich, y = virus_rich))+
  geom_point(aes(color = general_env_feature))+
  geom_smooth(method = "lm", aes(color = general_env_feature), se=F)+
  stat_cor(method = "pearson", aes(color = general_env_feature))+
  theme_bw() +
  xlab("Bacterial species richness")+
  ylab("vOTUS richness")+
  scale_color_manual(values = c("grey70", "grey30"))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_voturich_bactrich_final.pdf", device = "pdf", width = 5, height = 5 , units = "in")

# by vegetation type
ggplot(md_all, aes(x = bact_rich, y = virus_rich))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Bacterial richness")+
  ylab("vOTU richness")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_voturich_bactrich_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# PTUs vs vOTUs
ggplot(md_all, aes(x = plasmid_rich, y = virus_rich))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("Plasmid richness")+
  ylab("vOTUS richness")+
  facet_wrap(~general_env_feature, scales = "free")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_voturich_pturich.pdf", device = "pdf", width = 7, height = 5 , units = "in")

# by vegetation type
ggplot(md_all, aes(x = plasmid_rich, y = virus_rich))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("PTU richness")+
  ylab("vOTU richness")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_voturich_pturich_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# final corr plot
ggplot(md_all, aes(x = plasmid_rich, y = virus_rich))+
  geom_point(aes(color = general_env_feature))+
  geom_smooth(method = "lm", aes(color = general_env_feature), se=F)+
  stat_cor(method = "pearson", aes(color = general_env_feature))+
  theme_bw() +
  xlab("PTU richness")+
  ylab("vOTUS richness")+
  scale_color_manual(values = c("grey70", "grey30"))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_voturich_pturich_final.pdf", device = "pdf", width = 5, height = 5 , units = "in")





# ------------------------------------------------------------------------------------------------------------------------------
## correlations to ABIOTIC FACTORS

ggplot(md_all, aes(x = ags, y = exc_magnesium))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("vOTUs richness")+
  ylab("Magnesium")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_virusrich_mg_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


ggplot(md_all, aes(x = ags, y = exc_magnesium))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("PTU richness")+
  ylab("pH")+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')




# ------------------------------------------------------------------------------------------------------------------------------
# correlations to TRAITS
#md_traits <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_traits.csv")

#md_traits_tojoin <- md_traits %>%
  select(c("Sample", "ags", "acn", "gc_mean"))

#md_all_traits <- md_all %>%
  left_join(., md_traits_tojoin, by = "Sample")

  
###############
### LENGTH
############### 
# Plasmid length vs. AGS
ggplot(md_all, aes(x = av_len, y = ags))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("PTU size (bp)")+
  ylab("Average genome size (AGS)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ptulen_ags.pdf", device = "pdf", width = 7, height = 5 , units = "in")

# final corr plot
ggplot(md_all, aes(x =  gc_mean, y = ags))+
  geom_point(aes(color = general_env_feature))+
  geom_smooth(method = "lm", aes(color = general_env_feature), se=F)+
  stat_cor(method = "pearson", aes(color = general_env_feature))+
  theme_bw() +
  xlab("PTU size (bp)")+
  ylab("Average genome size (AGS)")+
  scale_color_manual(values = c("grey70", "grey30"))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ptulen_ag_final.pdf", device = "pdf", width = 5, height = 5 , units = "in")


# Plasmid length vs. AGS vegetation
ggplot(md_all, aes(x = av_len, y = ags))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("PTU size (bp)")+
  ylab("Average genome size (AGS)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ptulen_ags_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")

# Viral rich vs. AGS
ggplot(md_all, aes(x = virus_av_len, y = ags))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Viral genome size (bp)")+
  ylab("Average genome size (AGS)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_viruslen_ags_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")



# Plasmid length vs. proportion of broad range plasmids
ggplot(md_all, aes(x = gc_mean, y = virus_rich))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("PTU size (bp)")+
  ylab("Proportion of broad range plasmids")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ptulen_prop_broad.pdf", device = "pdf", width = 7, height = 5 , units = "in")



# Plasmid length vs. ACN
ggplot(md_all, aes(x = av_len, y = acn))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("PTU size (bp)")+
  ylab("Average copy number (ACN)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ptulen_acn_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")



# Viral length vs. ACN
ggplot(md_all, aes(x = virus_av_len, y = acn))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Viral genome size (bp)")+
  ylab("Average copy number (ACN)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_viruslen_acn_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")

ggplot(md_all, aes(x = virus_av_len, y = acn))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("Viral genome size (bp)")+
  ylab("Average copy number (ACN)")+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_viruslen_acn.pdf", device = "pdf", width = 7, height = 5 , units = "in")

# final corr plot
ggplot(md_all, aes(x =  virus_av_len, y = acn))+
  geom_point(aes(color = general_env_feature))+
  geom_smooth(method = "lm", aes(color = general_env_feature), se=F)+
  stat_cor(method = "pearson", aes(color = general_env_feature))+
  theme_bw() +
  xlab("Viral genome size (bp)")+
  ylab("Average copy number (ACN)")+
  scale_color_manual(values = c("grey70", "grey30"))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_viruslen_acn_final.pdf", device = "pdf", width = 5, height = 5 , units = "in")



# virus vs plasmid length
ggplot(md_all, aes(x = virus_av_len, y = av_len))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Viral genome size (bp)")+
  ylab("Plasmid size (bp)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'none')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_viruslen_ptlen_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


ggplot(md_all, aes(x = virus_av_len, y = av_len))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("Viral genome size (bp)")+
  ylab("Average genome size (AGS)")+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')


###############
### RELATIVE ABUNDANCE
############### 

# Plasmid relative abundance vs. AGS and ACN
ggplot(md_all, aes(x = plasmids_relab, y = acn))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Plasmid relative abundance")+
  ylab("Average copy number (ACN)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_plasmidrelab_acn_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# Virus relative abundance vs. AGS and ACN
ggplot(md_all, aes(x = virus_relab, y = ags))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type,fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Virus relative abundance")+
  ylab("Average genome size (AGS)")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_virusrelab_ags_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")

ggplot(md_all, aes(x = plasmids_relab, y = gc_mean))+
  geom_point()+
  geom_smooth(method = "lm")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("Virus relative abundance")+
  ylab("Average copy number (ACN)")+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')




### Proportion of broad range plasmids cs. ptu size

# final corr plot
ggplot(md_all, aes(x =  av_len, y = prop_broad))+
  geom_point(aes(color = general_env_feature))+
  geom_smooth(method = "lm", aes(color = general_env_feature), se=F)+
  stat_cor(method = "pearson", aes(color = general_env_feature))+
  theme_bw() +
  xlab("Plasmid size (bp)")+
  ylab("Proportion of broad range plasmids")+
  scale_color_manual(values = c("grey70", "grey30"))+
  theme(legend.position = "none", text = element_text(size = 16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_ptulen_propbroad_final.pdf", device = "pdf", width = 5, height = 5 , units = "in")




md_count <- md_all %>%
  group_by(vegetation_type) %>%
  count()






# ------------------------------------------------------------------------------------------------------------------------------
## Mantel tests

# To run this you have to load the distance matrices for each thing on their script

mantel_result <- mantel(bact.bray, plasmids.bray, method = "pearson", permutations = 999)
mantel_result <- mantel(bact.bray_arid, plasmids.bray_arid, method = "pearson", permutations = 999)
mantel_result <- mantel(bact.bray_temperate, plasmids.bray_temperate, method = "pearson", permutations = 999)

# Bacteria vs. PTUs
# extract pairwise differences
distances_df <- data.frame(
  Dist1 = as.vector(bact.bray),
  Dist2 = as.vector(plasmids.bray)
)

# Scatter plot with regression line
ggplot(distances_df, aes(x = Dist1, y = Dist2)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(
    title = paste("Mantel Test Plot (r =", round(mantel_result$statistic, 3), ", p =", mantel_result$signif, ")"),
    x = "Bacterial community composition",
    y = "Plasmid community composition"
  ) +
  theme_bw()
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/mantel_plasmid_bact_temperate.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# Bacteria vs. vOTUs
mantel_result <- mantel(bact.bray, virus.bray, method = "pearson", permutations = 999)
mantel_result <- mantel(bact.bray_arid, virus.bray_arid, method = "pearson", permutations = 999)
mantel_result <- mantel(bact.bray_temperate, virus.bray_temperate, method = "pearson", permutations = 999)

# extract pairwise differences
distances_df <- data.frame(
  Dist1 = as.vector(bact.bray),
  Dist2 = as.vector(virus.bray)
)

# Scatter plot with regression line
ggplot(distances_df, aes(x = Dist1, y = Dist2)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(
    title = paste("Mantel Test Plot (r =", round(mantel_result$statistic, 3), ", p =", mantel_result$signif, ")"),
    x = "Bacterial community composition",
    y = "Viral community composition"
  ) +
  theme_bw()
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/mantel_virus_bact_temperate.pdf", device = "pdf", width = 7, height = 5 , units = "in")


# PTUs vOTUs
mantel_result <- mantel(plasmids.bray, virus.bray, method = "pearson", permutations = 999)
mantel_result <- mantel(plasmids.bray_arid, virus.bray_arid, method = "pearson", permutations = 999)
mantel_result <- mantel(plasmids.bray_temperate, virus.bray_temperate, method = "pearson", permutations = 999)
# extract pairwise differences
distances_df <- data.frame(
  Dist1 = as.vector(plasmids.bray),
  Dist2 = as.vector(virus.bray)
)

# Scatter plot with regression line
ggplot(distances_df, aes(x = Dist1, y = Dist2)) +
  geom_point(color = "blue", alpha = 0.7) +
  geom_smooth(method = "lm", color = "red", se = FALSE) +
  labs(
    title = paste("Mantel Test Plot (r =", round(mantel_result$statistic, 3), ", p =", mantel_result$signif, ")"),
    x = "Plasmid community composition",
    y = "Virus community composition"
  ) +
  theme_bw()
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/mantel_plasmid_virus_arid.pdf", device = "pdf", width = 7, height = 5 , units = "in")




# --------------------------------------------------------------------------------------------------------------------
##### Corplot
md_arid <- md_all %>%
  filter(general_env_feature == "Arid")

md_temperate <- md_all %>%
  filter(general_env_feature == "Temperate")

# CORRPLOT
md_corplot <- md_arid %>% # change between arid and temperate
  dplyr::select(c("bact_rich", "plasmid_rich", "virus_rich", 
                  "ags", "acn", "gc_mean", "CRISPR_i",
                  "av_len","prop_broad", "virus_av_len", "ratio_vt"))



### Without bonferroni correction
testRes = corrplot::cor.mtest(md_corplot, conf.level = 0.95)
corrplot::corrplot(cor(md_corplot, use = "complete.obs"),
                   method = 'square',
                   type = "upper", 
                   sig.level = 0.05,           # Significance level
                   insig = "blank", 
                   tl.col = "black", 
                   p.mat = testRes$p,
                   pch.cex = 0.9, pch.col = 'grey20',
                   diag = FALSE)




# --------------------------------------------------------------------------------------------------------------------
##### Bonferroni correction


###################################################
#####  with p.adjust
# First, ensure we have a proper data frame with numeric values
md_fixed <- as.data.frame(lapply(md_corplot, unlist))

# Calculate the correlation matrix
cor_matrix <- cor(md_fixed, use="pairwise.complete.obs")

# Calculate raw p-values
p_values <- matrix(NA, ncol(md_fixed), ncol(md_fixed))
colnames(p_values) <- rownames(p_values) <- colnames(md_fixed)

for(i in 1:(ncol(md_fixed)-1)) {
  for(j in (i+1):ncol(md_fixed)) {
    test <- cor.test(md_fixed[,i], md_fixed[,j])
    p_values[i,j] <- test$p.value
    p_values[j,i] <- test$p.value  # Make it symmetric
  }
}

# Extract the p-values (excluding the diagonal)
p_values_vector <- p_values[upper.tri(p_values)]

# Apply Bonferroni correction using p.adjust
#adjusted_p_values <- p.adjust(p_values_vector, method = "bonferroni")
adjusted_p_values <- p.adjust(p_values_vector, method = "fdr")

# Put the adjusted p-values back into a matrix
adjusted_p_matrix <- matrix(NA, ncol(md_fixed), ncol(md_fixed))
colnames(adjusted_p_matrix) <- rownames(adjusted_p_matrix) <- colnames(md_fixed)
adjusted_p_matrix[upper.tri(adjusted_p_matrix)] <- adjusted_p_values

# Make the matrix symmetric
for(i in 1:(ncol(md_fixed)-1)) {
  for(j in (i+1):ncol(md_fixed)) {
    adjusted_p_matrix[j,i] <- adjusted_p_matrix[i,j]
  }
}

# Now create a correlation plot
library(corrplot)

corrplot(cor_matrix, 
         p.mat = adjusted_p_matrix,  # Use Bonferroni-adjusted p-values
         sig.level = 0.05,           # Significance level
         insig = "blank",            # Hide insignificant correlations
         method = "square",          # Visualization method
         type = "upper",             # Show only upper triangle
         tl.col = "black",           # Text label color
         diag = FALSE)               # Hide diagonal
