# Mobile OG analyses for base project


library(tidyverse)
library(vegan)
library(see)
library(ggdist)
library(ggplot2)
library(gghalves)  

pal1 <- c("orange2", "gold", "forestgreen")
pal1 <- c("#DC267F", "#FFB000", "#648FFF")
pal1 <- c("#D81B60", "#FFC107", "#004D40")

mges_counts <- read.table("/Volumes/BunnyBike/mge_urban/base/mobileOG/mges_count_table.txt", sep = "\t", header = T)
md <- read.table("/Volumes/BunnyBike/mge_urban/base/base_metadata.txt", sep = "\t", header = T)
mobileog_diamond <- read.table("/Volumes/BunnyBike/mge_urban/base/mobileOG/base_combined_mobileOG__121124.out")


## DIAMOND OUTPUT
colnames(mobileog_diamond) <- c("qseqid", "sseqid", "pident", "length", "mismatch", "gapopen", "qstart", 
                                "qend", "sstart", "send", "evalue", "bitscore")
# keep only the best hit per gene
mobileog_diamond_filtered <- mobileog_diamond %>%
  dplyr::group_by(qseqid) %>% # group by same read name
  dplyr::arrange(bitscore,.by_group = T) %>% # sort by evalue within group
  slice_head(n = 1) %>% # select the first hit (best one) for each group 
  filter(!pident < 90) %>%
  filter(length > 50)

## TAKES VERY LONG! 

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
  filter(!MGE_type == "Plasmid") %>%
  mutate(Contig = gsub(qseqid, pattern = "_[0-9]*$", replacement = ""))


## COUNTS OUTPUT

# keep only the mges that we filtered out from the diamond output
counts_nonplasmids <- mges_counts %>%
  filter(Contig %in% non_plasmid_diamond$Contig) 
# no change from MGEs counts

#####
# ADD MOBILE_OG INFORMATION
#####

## select blast information we want
non_plasmid_diamont_tojoin <- non_plasmid_diamond %>%
  mutate(Contig = gsub(qseqid, pattern = "_[0-9]$", replacement = "")) %>%
  select(c("Contig", "MGE_type", "mobileOG_id"))


counts_nonplasmids_final <- counts_nonplasmids %>%
  left_join(.,non_plasmid_diamont_tojoin, by = "Contig")

counts_nonplasmids_final_grouped <- counts_nonplasmids %>%
  left_join(.,non_plasmid_diamont_tojoin, by = "Contig") %>%
  # GROUP BY MOBILEOG_ID
  select(-c("Contig", "MGE_type", "qseqid")) %>%
  group_by(mobileOG_id) %>%
  summarize_all(sum)



mges_counts_long <- counts_nonplasmids_final_grouped %>%.  ### NOTICE HERE I AM USING THE GROUPED MGES! 
  pivot_longer(-mobileOG_id, names_to = "Sample", values_to = "Counts") %>%
  mutate(Sample = gsub(Sample, pattern = ".Read.Count", replacement = "")) %>%
  mutate(Sample = gsub(Sample, pattern = "X", replacement = ""))

md_filtered <- md %>%
  mutate(Sample = gsub(sample_id, pattern = "_[0-9]", replacement = "")) %>%
  filter(Sample %in% mges_counts_long$Sample) %>%
  filter(!general_env_feature == "Australia") %>%
  filter(!is.na(vegetation_type)) %>%
  filter(!vegetation_type == "Forest")  

mges_counts_filtered <- mges_counts_long %>%
  filter(Sample %in% md_filtered$Sample) %>%
  pivot_wider(names_from = mobileOG_id, values_from = Counts ) %>%
  arrange(desc(Sample)) %>%
  column_to_rownames(var = "Sample") 


md_final<- md_filtered %>%
  arrange(desc(Sample)) %>%   # check that the
  column_to_rownames(var = "Sample") 




# make sure the rows are in the same order
rownames(mges_counts_filtered) == rownames(md_final)

summary(rowSums(mges_counts_filtered))

# richness rarified
md_final$mges_rich <- specnumber(rrarefy(mges_counts_filtered, sample = 1960))


ggplot(md_final, aes(x = vegetation_type, y = mges_rich))+
  #geom_violin(aes(fill =urban.natural, color = urban.natural), binwidth = 1, alpha=0.4)+
  geom_boxplot(alpha=0.6, outlier.shape = NA )+
  geom_jitter(aes(color= vegetation_type), size = 2.5)+
  xlab(NULL) + 
  ylab("MGE richness")+
  scale_fill_see()+
  scale_color_see()+
  theme_bw()+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))+
  facet_wrap(~general_env_feature)
ggsave("/Volumes/BunnyBike/mge/base/figures/mge_rich.pdf", device = "pdf", width = 5, height = 5 , units = "in")



# nice plot
ggplot(md_final, aes(x = vegetation_type, y = mges_rich, color = vegetation_type, fill = vegetation_type))+
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
  ylab("MGE richness")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/mge_rich_nice_grouped.pdf", device = "pdf", width = 7, height = 5 , units = "in")

anova(lm(mges_rich ~  vegetation_type * general_env_feature, data = md_final))

# richness as number of mge contigs/number of reads
md_final <- md_final %>%
  mutate(mge_rich_2 = rowSums(mges_counts_filtered)/tot_reads)


ggplot(md_final, aes(x = vegetation_type, y = mge_rich_2, color = vegetation_type, fill = vegetation_type))+
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
  ylab("MGE diversity")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/mge_rich_2_nice_grouped.pdf", device = "pdf", width = 7, height = 5 , units = "in")


anova(lm(mge_rich_2 ~  vegetation_type * general_env_feature, data = md_final))

## ## ## ## ## ## ## ## ## 
## RPKM
## ## ## ## ## ## ## ## ## 
# 1. prepare the length and the total reads to left join
# import and format the table of lengths per gene
mges_len <- read_csv("/Volumes/BunnyBike/mge_urban/local/mobileOG/mobileOG-db_beatrix-1.6.All_length.csv")
colnames(mges_len) <- c("seq_name", "len")

mges_len_fixed <- mges_len %>%
  separate(seq_name, into = c("mobileOG_id", "gene_name", "other_id", "category", "something", "db", "class"), sep = "\\|") %>%
  mutate(len = gsub(len, pattern = "[A-Za-z]", replacement = "")) %>%
  mutate(len = gsub(len, pattern = "\\|", replacement = "")) %>%
  mutate(len = gsub(len, pattern = ",", replacement = "")) %>%
  mutate(len = gsub(len, pattern = " ", replacement = "")) %>%
  mutate(len = gsub(len, pattern = "\\/", replacement = "")) %>%
  mutate(len = gsub(len, pattern = "\\-", replacement = "")) %>%
  select(c("mobileOG_id", "len")) %>%
  mutate(len_num = as.numeric(len)) %>%
  select(-c("len"))

# select the total reads per sample from the metadata
total_reads <- md_final %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(sample, tot_reads)

# 2. pivot longer the count table to make every gene in every sample a row
# USING MOBILEOG_ID
mges_long <- mges_counts_filtered %>%
  rownames_to_column(var = "sample") %>%
  pivot_longer(cols = -c("sample"), values_to = "Count", names_to= "mobileOG_id") %>%
  # group the genes with the same mobileOG_id and sum their counts
  dplyr::group_by(mobileOG_id, sample) %>%
  dplyr::summarize(Counts = sum(Count))


# 3. join the total reads to each row based on sample
mges_long_plusreads <- left_join(mges_long, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
mges_long_plusreads_plustlen <- left_join(mges_long_plusreads, mges_len_fixed, by = "mobileOG_id")
#select(-c("mobileOG_id")) %>% # for when using gene
#unique()


# 5. compute the RPKM calculation and create count table from it 
# USING MOBILEOG_ID
rpkm_mges_final <- mges_long_plusreads_plustlen  %>%
  dplyr::mutate(num = (Counts*10e6)) %>%
  dplyr::mutate(denom = tot_reads*len_num) %>%
  dplyr::mutate(RPKM = num/denom) %>%
  dplyr::select(-c("Counts", "tot_reads", "len_num", "num", "denom")) %>%
  pivot_wider(names_from = sample, values_from = RPKM) %>%
  column_to_rownames(var = "mobileOG_id") %>%
  drop_na()



'''
# RPKM USING GENES
mges_long_gene <- counts_nonplasmids_final %>%
  select(-c("MGE_type")) %>%
  pivot_longer(cols = -c("Gene", "mobileOG_id"), values_to = "Count", names_to= "sample") #%>%

# 3. join the total reads to each row based on sample
mges_long_plusreads_gene <- left_join(mges_long_gene, total_reads, by = "sample")
# 4. join the gene length to each row based on gene
mges_long_plusreads_plustlen_gene <- left_join(mges_long_plusreads_gene, mges_len_fixed, by = "mobileOG_id")
mges_long_plusreads_plustlen_gene <- mges_long_plusreads_plustlen_gene %>%
  filter(!sample == "Ex45-T1-9") #%>%
#select(-c("mobileOG_id")) %>% # for when using gene
#unique()

# 5. compute the RPKM calculation and create count table from it 
# USING GENES
rpkm_mges_final_gene <- mges_long_plusreads_plustlen_gene  %>%
  dplyr::mutate(RPKM = (as.numeric(Count)*10e6)/(as.numeric(total.reads)*as.numeric(len))) %>%
  dplyr::select(-c("Count", "total.reads", "len")) %>%
  #group_by(mobileOG_id, sample) %>%
  #dplyr::summarize(mean_RPKM = mean(RPKM), ) %>%
  pivot_wider(names_from = sample, values_from = RPKM) %>%
  column_to_rownames(var = "Gene") %>%
  drop_na()

'''


############################
# abundance
############################
sum_rel_ab_mge <- as.data.frame(t(rpkm_mges_final)) 

sum_rel_ab_mge_final <- sum_rel_ab_mge %>%
  #mutate_all(., function(x) as.numeric(as.character(x))) %>%
  dplyr::mutate(mge_relab = rowSums(across(where(is.numeric))))%>%
  dplyr::select(mge_relab) %>%
  rownames_to_column(var = "sample")

md_final <- md_final %>%
  rownames_to_column(var = "sample") %>%
  left_join(., sum_rel_ab_mge_final, by = "sample")

# nice plot
ggplot(md_final, aes(x = vegetation_type, y = mge_relab, color = vegetation_type, fill = vegetation_type))+
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
  ylab("MGE abundance RPKM")+
  theme_bw()+
  facet_wrap(~general_env_feature)+
  theme(legend.position = "none", text = element_text(size=16), axis.text.x = element_text(angle = 45, hjust = 1))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/mge_relab_nice_grouped.pdf", device = "pdf", width = 7, height = 5 , units = "in")

anova(lm(mge_relab~  vegetation_type * general_env_feature, data = md_final))

############################
# beta-diversity
############################
mges.bray <- vegdist(t(rpkm_mges_final), method="bray")

#ordination (non-multidimensional scaling)
mges.nmds <- metaMDS(mges.bray, k=2, try = 100)
md_final$Axis01 = mges.nmds$points[,1]
md_final$Axis02 = mges.nmds$points[,2]
mges.nmds$stress #0.03989818 the smaller the better, good <0.3)

ggplot(md_final, aes(Axis01, Axis02))+
  geom_point(aes(color=vegetation_type, shape = general_env_feature), size=4)+
  #stat_ellipse(aes(color = vegetation_type))+
  stat_ellipse(aes(group = general_env_feature, linetype = general_env_feature))+
  scale_shape_manual(values=c(16, 13))+
  scale_color_manual(values = pal1)+
  theme_bw()+
  theme(legend.position="bottom", text = element_text(size=12))
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/mges_beta_rpkm_grouped.pdf", device = "pdf", width = 7, height = 5 , units = "in")

adonis2(mges.bray ~ vegetation_type * general_env_feature, data = md_final,  permutations = 999, method = "bray")




