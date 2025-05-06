# plasmid taxonomy BASE project

library(tidyverse)
library(vegan)


hotspot <- read_tsv("/Volumes/BunnyBike/mge_urban/base/genomad_all_contigs/host_lineage.tsv")



# barplot of phyhla

phyla <- hotspot %>%
  group_by(Phylum) %>%
  summarize(count = n()) %>%
  mutate(pepito = "pepito") %>%
  mutate(count = as.numeric(count)) %>% 
  arrange(desc(count))

# barplot

pal3 <- c("grey",  "#FEE08B", "#3288BD", "#D53E4F" ,"#ABDDA4", "purple","#9E0142","#FDAE61" , "#F46D43", "pink", "black", "darkgreen")

# reorder

ggplot(phyla, aes(x = pepito, y = count, fill = Class, group = count))+
  geom_col()+
  xlab("")+
  ylab("Number of plasmids")+
  scale_fill_manual(values = pal3)+
  theme_bw()+
  theme(axis.text.x = element_blank())
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/host_taxonomy.pdf", device = "pdf", width = 4, height = 5 , units = "in")



# PIE PLOT
phyla_pie <- phyla %>%
  group_by(Phylum) %>%
  summarise(
    total_count = sum(count),
    percentage = round(total_count / sum(phyla$count) * 100, 2)
  ) %>%
  mutate(
    # Add percentage sign
    label = paste0(round(percentage, 1), "%"),
    # Calculate the position for labels
    position = cumsum(percentage) - percentage/2, 
    pos_rad = position * pi / 180
  )

library(see)
ggplot(phyla_pie , aes(x = 2, y = percentage, fill = Phylum)) +
  # Create the donut
  geom_bar(stat = "identity", width = 1) +
  # Add percentage labels
  geom_text(aes(label = label),
            position = position_stack(vjust = 0.5)) +
  geom_text(aes(label = Phylum),
            position = position_stack(vjust = 0.4)) +
  # Convert to polar coordinates to make it circular
  coord_polar(theta = "y", start = 0) +
  # Remove unnecessary elements and customize theme
  theme_void() +
  #scale_fill_see()+
  # Create a hole in the center to make it a donut
  xlim(0.5, 2.7) +
  # Custom colors and legend position
  theme(legend.position = "bottom", plot.margin = margin(1, 1, 1, 1, "cm"))



# Alluvial with claude
# First, let's load required libraries
library(ggplot2)
library(tidyr)
library(dplyr)
library(ggtree)
library(ggalluvial)

# Create the alluvial plot with a modified approach
taxonomy_alluvial <- hotspot %>%
  # First, store the Phylum information before gathering
  mutate(Original_Phylum = Phylum) %>%
  gather(key = "Level", value = "Taxon", Phylum:Species) %>%
  mutate(
    Level = factor(Level, levels = c("Phylum", "Class", "Order", "Family", "Genus", "Species")),
    Taxon = ifelse(Taxon == "-", "Unclassified", Taxon)
  )

# Create the plot
ggplot(taxonomy_alluvial,
       aes(x = Level, stratum = Taxon, alluvium = Contig, fill = Original_Phylum)) +
  geom_alluvium(width = 1/12) +
  geom_stratum(width = 1/12, fill = "white", color = "grey") +
  theme_minimal() +
  geom_text(stat = "stratum", aes(label = Taxon), size = 3, hjust = 1) +  # Add labels
  theme() +
  labs(
    y = "Number of PTUs",
    fill = "Phylum") +
  scale_fill_brewer(palette = "Accent")  # Optional: adds a nice color palette



# ------------------------------------------------------------------------------------------------------------
# Plot taxonomy hosts by abundance plasmids each sample
# NOTE: run plasmids_prelim_analyses.R first to get md and rpkm data

md_join <- md_final_final %>%
  rownames_to_column(var = "sample") %>%
  dplyr::select(c("sample", "general_env_feature", "vegetation_type"))

# Join the pab table to the taxonomy
taxonomy_pab <- rpkm_plasmids_final %>%
  rownames_to_column(var = "Contig") %>%
  left_join(., hotspot, by = "Contig") %>%
  dplyr::select(-c("Contig", "Class", "Order", "Family", "Genus", "Species")) %>%
  pivot_longer(-Phylum, names_to = "sample", values_to = "pab") %>%
  left_join(., md_join, by = "sample") %>%
  group_by(Phylum, general_env_feature, vegetation_type) %>%
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

pal3 <- c("#FEE08B" ,"olivedrab", "#D53E4F" ,"orange", 
          "purple","chocolate4","#3288BD", "navy")

ggplot(taxonomy_pab, aes(x = vegetation_type, y = percentage, fill = Phylum, group = counts))+
  geom_col()+
  scale_fill_manual(values = pal3)+
  xlab("")+
  ylab("Proportion of counts (%)")+
  facet_wrap(~general_env_feature)+
  theme_bw()+
  theme(legend.position = "none")
ggsave("/Volumes/BunnyBike/mge_urban/base/figures/plasmid_hosts_phyla_nolegend.pdf", device = "pdf", width = 6, height = 5 , units = "in")


ggplot(taxonomy_pab, aes(x = vegetation_type, y = percentage, fill = Phylum, group = counts))+
  geom_col()+
  scale_fill_manual(values = pal3)+
  xlab("")+
  ylab("Proportion of counts (%)")+
  facet_wrap(~general_env_feature)+
  theme_bw()+
  theme(legend.position = "none")
