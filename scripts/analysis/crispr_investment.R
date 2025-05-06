# CRISPR investment calculation and visualization for base project


library(tidyverse)
library(car)
library(lme4)

pal1 <- c("#D81B60", "#997404", "#004D40")

# load crispr results
crispr_results <- read_tsv("/Volumes/BunnyBike/mge_urban/base/traits/crispr_spacer_counts.tsv") %>%
  mutate(Sample = as.numeric(Sample))

# already has genome number in it
md_all <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_traits.csv")

# join and calculte crispr investment
md_final <- md_all %>%
  left_join(., crispr_results, by = "Sample") %>%
  mutate(CRISPR_i = Num_Spacers/genome_number)
  

md_final$combined_var <- paste(md_final$general_env_feature, md_final$vegetation_type)
md_final$combined_var <- factor(md_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                  "Arid Shrubland", "Temperate Shrubland", 
                                                                  "Arid Woodland", "Temperate Woodland"))

ggplot(md_final, aes(x = vegetation_type, y = CRISPR_i, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  ylab("CRISPR investment")+ 
  theme_bw()+
  #facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/crispr_investment_paralell.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(CRISPR_i ~  vegetation_type * general_env_feature + (1|season), data = md_final), type = 3)
#vegetation_type                     0.3555  2   0.837138   
#general_env_feature                 2.3696  1   0.123717   
#vegetation_type:general_env_feature 1.5202  2   0.467608 
r2(lmer(CRISPR_i ~  vegetation_type * general_env_feature + (1|season), data = md_final))

write_csv(md_final, "/Volumes/BunnyBike/mge_urban/base/md_final_crispr.csv")

### Correlations
ggplot(md_final, aes(x = av_len, y = CRISPR_i))+
  geom_point()+
  geom_smooth(method = "lm", color = "black")+
  stat_cor(method = "pearson")+
  theme_bw() +
  xlab("Plasmid average length")+
  ylab("CRISPR_i")+
  facet_wrap(~general_env_feature, scales = "free")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_plasmidlength_crispr.pdf", device = "pdf", width = 7, height = 5 , units = "in")


ggplot(md_final, aes(x = av_len, y = CRISPR_i))+
  geom_point(aes(color = vegetation_type))+
  geom_smooth(method = "lm", aes(color = vegetation_type, fill = vegetation_type))+
  stat_cor(method = "pearson", aes(color = vegetation_type))+
  theme_bw() +
  xlab("Plasmid average length")+
  ylab("CRISPR_i")+
  scale_color_manual(values=pal1)+
  scale_fill_manual(values=pal1)+
  facet_wrap(~general_env_feature, scales = "free")+
  theme(legend.position = 'bottom')
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/cor_plasmidlength_crispr_veg.pdf", device = "pdf", width = 7, height = 5 , units = "in")




md_corplot <- md_final %>%
  dplyr::select(c("bact_rich", "plasmid_rich", "virus_rich", 
                  "ags", "acn", "gc_mean", "CRISPR_i", 
                  "av_len", "virus_av_len", 
                  "ph", "organic_carbon", "exc_magnesium"))

testRes = corrplot::cor.mtest(md_corplot, conf.level = 0.95)
corrplot::corrplot(cor(md_corplot, use = "complete.obs"),
                   method = 'square',
                   type = "upper", 
                   tl.col = "black", 
                   addCoef.col = 'black',
                   p.mat = testRes$p,
                   insig =  "label_sig",
                   sig.level = c(0.001, 0.01, 0.05),
                   pch.cex = 0.9, pch.col = 'grey20',
                   diag = FALSE)






