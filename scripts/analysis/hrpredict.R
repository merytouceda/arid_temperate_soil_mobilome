# analysis of the Host range data of BASE

library(tidyverse)
library(vegan)
pal1 <- c("#D81B60", "#FFC107", "#004D40")
pal1 <- c("#D81B60", "#997404", "#004D40")

# Load and clean 
host_range_s <- read_tsv("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/HostRange_Species.tsv")

separations <- as.character(c(1:205))

host_range_s_clean <- host_range_s %>%
  add_row(!!!setNames(colnames(.), colnames(.)), .before = 1) %>%
  # Give new generic column names
  set_names(c("plasmid", paste0("col", 2:ncol(.)))) %>%
  # Remove any asterisks from the first row values
  mutate(across(everything(), ~if_else(row_number() == 1, 
                                       str_remove_all(., "\\*"), 
                                       .))) %>%
  mutate(col2 = gsub(col2, pattern = "f__", replacement = "\\| f__")) %>%
  mutate(col2= gsub(col2, pattern = "^\\|", replacement = "")) %>%
  filter(!col2 == "No species was predicted") %>%
  mutate(degree = case_when(str_detect(col3, pattern = "f__") ~ 'IV, V, VI', 
                            TRUE ~ 'I, II, III'))



##### PLOT PRESENCE OF DEGREE IN URBAN AND NATURAL

# Join this information with presence absence in urban and natural

ptu_counts <- read.table("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/combined_assembly/ptus_count_table.txt", sep = "\t", header = T)
md <- read_csv("/Volumes/BunnyBike/mge_urban/base/md_final_final_traits.csv")

plasmids_counts_long <- ptu_counts  %>%
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
  arrange(desc(Sample)) %>%   # check that the
  column_to_rownames(var = "Sample") 


md_final$combined_var <- paste(md_final$general_env_feature, md_final$vegetation_type)
md_final$combined_var <- factor(md_final$combined_var, levels = c("Arid Grassland", "Temperate Grassland",
                                                                  "Arid Shrubland", "Temperate Shrubland", 
                                                                  "Arid Woodland", "Temperate Woodland"))

md_to_join <- md_final %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(c("sample", "general_env_feature", "vegetation_type", "season", "combined_var"))


ptu_pab <- plasmids_counts_filtered %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(-c("sample"), names_to = "plasmid", values_to = "count") %>%
  mutate(sample = gsub(sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(sample = gsub(sample, pattern = "X", replacement = "")) %>%
  mutate(pab = case_when(count ==  0 ~ 0, 
                         TRUE ~ 1)) %>%
  dplyr::select(-c("count")) %>%
  pivot_wider(names_from = "sample", values_from = "pab")

host_range_count <- host_range_s_clean %>%
  dplyr::select(c("plasmid","degree")) %>%
  left_join(., ptu_pab, by = "plasmid") %>%
  drop_na() %>%
  dplyr::select(-plasmid) %>%
  pivot_longer(-degree, names_to = "sample", values_to="count") %>%
  group_by(degree, sample) %>%
  summarize(counts = sum(count)) %>%
  left_join(., md_to_join, by = "sample")



ggplot(host_range_count, aes(x = vegetation_type, y = counts, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  ylab("Number of plasmids")+ 
  theme_bw()+
  facet_wrap(~degree)+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/host_range_parallel.pdf", device = "pdf", width = 7, height = 5 , units = "in")

# STATS??
host_range_count_broad <- host_range_count %>%
  filter(degree == "IV, V, VI")

Anova(lmer(counts ~  vegetation_type * general_env_feature  + (1|season), data = host_range_count_broad ), type = 3)

# calculate percentage
degree_veg_percentage <- host_range_count %>%
  group_by(degree, general_env_feature) %>%
  summarise(total_counts = sum(counts), .groups = "drop") %>%
  group_by(general_env_feature) %>%
  mutate(degree_total = sum(total_counts),
         within_degree_percentage = total_counts / degree_total * 100) %>%
  ungroup() %>%
  mutate(overall_percentage = total_counts / sum(total_counts) * 100) %>%
  mutate(within_degree_percentage = round(within_degree_percentage, 2),
         overall_percentage = round(overall_percentage, 2))


# separate arid and temperate
host_range_count_arid <- host_range_count %>%
  filter(general_env_feature == "Arid")

host_range_count_temp<- host_range_count %>%
  filter(!general_env_feature == "Arid")


ggplot(host_range_count_arid, aes(x = vegetation_type, y = counts))+
  geom_boxplot(aes(color = vegetation_type))+
  geom_jitter(aes(color = vegetation_type))+
  facet_wrap(~degree)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Number plasmids")+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/host_range_arid.pdf", device = "pdf", width = 8, height = 8 , units = "in")

ggplot(host_range_count_temp, aes(x = vegetation_type, y = counts))+
  geom_boxplot(aes(color = vegetation_type))+
  geom_jitter(aes(color = vegetation_type))+
  facet_wrap(~degree)+
  scale_fill_manual(values=pal1)+
  scale_color_manual(values=pal1)+
  xlab(NULL) + 
  ylab("Number plasmids")+
  theme_classic()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/host_range_temp.pdf", device = "pdf", width = 8, height = 8 , units = "in")


# calculate proportion of broad host range plasmids

host_range_ratio <- host_range_count %>%
  ungroup() %>%
  select(c("sample", "degree", "counts")) %>%
  mutate(degree_term = case_when(degree == "IV, V, VI"~ "broad", 
                                 TRUE ~ "narrow")) %>%
  select(-degree) %>%
  pivot_wider(names_from = "degree_term", values_from = "counts") %>%
  mutate(prop_broad = broad/(narrow + broad)) %>%
  left_join(., md_to_join, by = "sample")

ggplot(host_range_ratio, aes(x = vegetation_type, y = prop_broad, color = combined_var, fill = combined_var))+
  geom_boxplot(width = .4, outlier.shape = NA, linewidth = 0.6, alpha = 0.3) +
  #ggdist::stat_halfeye(adjust = .33,  bandwidthwidth = .67, color = NA,  position = position_nudge(x = .15)) +
  #geom_point(position = position_nudge(x = .15))+
  gghalves::geom_half_point(side = "l", range_scale = .3, alpha = 1, size = 2.5)+
  #scale_fill_manual(values = pal1)+
  scale_fill_manual(values = c("lightpink", "#D81B60","#FFC107", "#997404", "green4","#004D40"))+
  #scale_color_manual(values=pal1)+
  scale_color_manual(values = c("lightpink", "#D81B60","#FFC107","#997404", "green4","#004D40"))+
  xlab(NULL) + 
  ylab("Proportion of broad range plasmids")+ 
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/proportion_broad_parallel.pdf", device = "pdf", width = 5, height = 5 , units = "in")

Anova(lmer(prop_broad ~  vegetation_type * general_env_feature  + (1|season), data = host_range_ratio), type = 3)
#vegetation_type                         1.8997  2    0.38680    
#general_env_feature                     3.3383  1    0.06768 .  
#vegetation_type:general_env_feature     1.5280  2    0.46580 
r2(lmer(prop_broad ~  vegetation_type * general_env_feature  + (1|season), data = host_range_ratio))

# calculate mean of broad range proportion
mean(host_range_ratio$prop_broad)
sd(host_range_ratio$prop_broad)

host_range_ratio_join <- host_range_ratio %>%
  select(c("sample", "prop_broad"))

md_final <- md_final %>%
  rownames_to_column(var = "sample") %>%
  left_join(host_range_ratio_join, by = "sample")

# save metadata for correlations
write_csv(md_final, "/Volumes/BunnyBike/mge_urban/base/md_final_hostrange.csv")

